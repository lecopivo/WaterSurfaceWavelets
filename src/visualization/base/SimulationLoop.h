#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

template <typename Simulation> class SimulationLoop {
public:
  ~SimulationLoop() {
    keep_running = false;
    main_loop_thread.join();
  }

  using Task = std::function<void(Simulation &)>;

  void start(Simulation &sim) {

    addTaskPerLoop([](Simulation &s) { s.update(); });

    // start the main thread
    main_loop_thread = std::thread([&sim, this]() {

      while (keep_running) {

        std::queue<Task> tasks;

        // Tasks per loop
        tasks = tasks_per_loop.getCopy();
        while (!tasks.empty()) {
          auto &task = tasks.front();
          task(sim);
          tasks.pop();
        }

        // Tasks to evaluate
        tasks = tasks_queue.getCopy();
        tasks_queue.clear();
        while (!tasks.empty()) {
          auto &task = tasks.front();
          task(sim);
          tasks.pop();
        }
      }
    });
  }

  void addTaskPerLoop(Task const &task) { tasks_per_loop.push(task); }

  void addTask(Task const &task) { tasks_queue.push(task); }

private:
  template <class T> class synch_queue {
  public:
    std::queue<T> getCopy() {
      std::lock_guard<std::mutex> guard(_mutex);
      return _queue;
    }

    void push(T const &t) {
      std::lock_guard<std::mutex> guard(_mutex);
      _queue.push(t);
    }

    void push(T &&t) {
      std::lock_guard<std::mutex> guard(_mutex);
      _queue.push(std::forward<T>(t));
    }

    void clear() {
      std::lock_guard<std::mutex> guard(_mutex);
      while (!_queue.empty())
        _queue.pop();
    }

  private:
    std::mutex    _mutex;
    std::queue<T> _queue;
  };

private:
  std::thread       main_loop_thread;
  synch_queue<Task> tasks_queue;
  synch_queue<Task> tasks_per_loop;

  bool keep_running = true;
};
