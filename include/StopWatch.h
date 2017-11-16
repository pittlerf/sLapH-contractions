#pragma once

/**
   Usage example:

   StopWatch sw;
   
   #pragma omp parallel
   {
       sw.start();
   
       // Heavy work
   
       sw.stop();
   }
   
   double wtime = sw.mean();
*/
class StopWatch {
  StopWatch(int const threads = -1)
      : starts_{0, threads == -1 ? omp_get_num_threads() : threads},
        ends_{0, starts_.size()} {}

  void start() { store_to(starts_); }

  void stop() { store_to(ends_); }

  double mean() const {
#pragma omp barrier
    double sum = 0.0;
    for (int i = 0; i < starts_.size(); ++i) {
      sum += ends_[i] - starts_[i];
    }
    return sum / starts_.size();
  }

 private:
  void store_to(std::vector<double> &vec) {
    int const t_id = omp_get_thread_num();
    double const wtime = omp_get_wtime();
    vec[t_id] = wtime;
  }

  std::vector<double> starts_;
  std::vector<double> ends_;
  bool started_ = false;
};

