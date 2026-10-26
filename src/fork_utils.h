#ifndef FLIPFISH_FORK_UTILS_H
#define FLIPFISH_FORK_UTILS_H

#include <Rcpp.h>
#include <vector>
#include <functional>
#include <string>

#ifndef _WIN32
#include <array>
#include <cerrno>
#include <cstring>
#include <unistd.h>    // fork, pipe, read, write, close
#include <sys/wait.h>  // waitpid
#endif

/*
 * Shared fork/pipe map-reduce helper for the handful of main.cpp functions that distribute
 * per-barcode work across processes (expected_bc_counts(), make_SpotSim(), make_SpotSim_DG()).
 * Replaces three near-identical copies of pipe setup / fork / write-loop / read-loop / waitpid
 * with one template, parameterized only by the flat payload element type T (double or int here)
 * and three caller-supplied callbacks:
 *
 *   batch_payload_size(i) -- size (element count) of batch i's flat result buffer, computed from
 *                            batch metadata alone (no work done), so the parent can size its read
 *                            buffer before the child's result arrives.
 *   compute_batch(i)      -- does batch i's work and returns its flat result buffer. Runs inside
 *                            the forked child for batch i (POSIX, n_forks > 1), or directly in the
 *                            calling process otherwise.
 *   reduce_batch(i, buf)  -- consumes batch i's buffer in the calling process (accumulate into the
 *                            caller's own result containers). Always runs in the parent/caller.
 *   batch_error_context(i) -- optional extra detail appended to the "unexpected EOF" error message
 *                            (e.g. which barcodes were in the batch); pass nullptr to omit.
 *
 * Windows has no fork(): MinGW/Rtools provides neither a working fork() nor pipe()/waitpid() with
 * POSIX semantics, so that code cannot even compile there. R's check_forks() (R/main.R) already
 * forces n_forks <- 1 at the call site on Windows, but that's a runtime guard and does nothing for
 * compilation -- the fork()/pipe() calls below still need to not exist in the Windows translation
 * unit at all. Hence the #ifdef _WIN32 around the whole parallel path: on Windows every batch just
 * runs serially in-process, regardless of n_forks, which is both correct (same result, no
 * parallelism) and sufficient, since the R layer never asks for more than one batch there anyway.
 */
template <typename T>
void run_forked_batches(
    int                                              n_forks,
    std::function<size_t(int)>                       batch_payload_size,
    std::function<std::vector<T>(int)>               compute_batch,
    std::function<void(int, const std::vector<T>&)>  reduce_batch,
    std::function<std::string(int)>                  batch_error_context = nullptr
  ) {

#ifdef _WIN32
  (void)batch_payload_size; // unused on this path: compute_batch's own result size is used directly
  (void)batch_error_context;
  for (int i = 0; i < n_forks; ++i) {reduce_batch(i, compute_batch(i));}
#else
  if (n_forks <= 1) {
    for (int i = 0; i < n_forks; ++i) {reduce_batch(i, compute_batch(i));}
    return;
  }

  // Pipes for inter-process communication
  std::vector<int>                pids(n_forks, -1);
  std::vector<std::array<int, 2>> pipes(n_forks);

  // Initialize pipes
  for (int i = 0; i < n_forks; ++i) {
    if (pipe(pipes[i].data()) != 0) {Rcpp::stop("Failed to create pipe!");}
  }

  // Fork processes
  for (int i = 0; i < n_forks; i++) {
    pid_t pid = fork();

    if (pid == 0) { // child process

      // Close unrelated pipe ends
      for (int j = 0; j < n_forks; ++j) {
        if (j == i) {
          close(pipes[j][0]); // keep write end
        } else {
          close(pipes[j][0]);
          close(pipes[j][1]);
        }
      }

      std::vector<T> buf = compute_batch(i);

      // Send result
      const char* buffer        = reinterpret_cast<const char*>(buf.data());
      size_t      nbytes        = sizeof(T) * buf.size();
      size_t      total_written = 0;
      while (total_written < nbytes) {
        ssize_t n_written = write(pipes[i][1], buffer + total_written, nbytes - total_written);
        if (n_written <= 0) {
          if (errno == EINTR) continue; // Retry if interrupted
          close(pipes[i][1]);
          _exit(1);
        }
        total_written += static_cast<size_t>(n_written);
      }

      close(pipes[i][1]); // Close write end
      _exit(0);            // Exit child process

    } else if (pid > 0) { // parent process
      pids[i] = pid;       // Grab child pid
      close(pipes[i][1]);  // Close write end
    } else {
      Rcpp::stop("Fork failed!");
    }
  }

  // Fetch results from pipes
  for (int i = 0; i < n_forks; i++) {
    size_t         payload_size = batch_payload_size(i);
    std::vector<T> buf(payload_size);

    // Read the row from the pipe into the buffer
    char*  buffer     = reinterpret_cast<char*>(buf.data());
    size_t nbytes     = sizeof(T) * payload_size;
    size_t total_read = 0;
    while (total_read < nbytes) {
      ssize_t n_read = read(pipes[i][0], buffer + total_read, nbytes - total_read);
      if (n_read == 0) {
        close(pipes[i][0]);
        int status;
        waitpid(pids[i], &status, WNOHANG);
        std::string msg =
          "\nUnexpected EOF while reading pipe " + std::to_string(i) + " of " + std::to_string(n_forks) + ". " +
          "\nCheck for quiet core crash from a child-process error." +
          "\nRead " + std::to_string(total_read) + " of " + std::to_string(nbytes) + " bytes.";
        if (batch_error_context) {msg += batch_error_context(i);}
        Rcpp::stop(msg);
      }
      if (n_read < 0) {
        if (errno == EINTR) continue; // Retry if interrupted
        close(pipes[i][0]);
        Rcpp::stop("Pipe read failed");
      }
      total_read += static_cast<size_t>(n_read);
    }

    reduce_batch(i, buf);

    close(pipes[i][0]); // Close read end
    int status;
    waitpid(pids[i], &status, 0);
    if (WIFSIGNALED(status)) {Rcpp::Rcout << "Child killed by signal " << WTERMSIG(status) << std::endl;}
    if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {Rcpp::Rcout << "Child exited with code " << WEXITSTATUS(status) << std::endl;}
  }
#endif
}

#endif // FLIPFISH_FORK_UTILS_H
