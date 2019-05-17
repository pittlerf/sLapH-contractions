#include <timings.hpp>

void baz(std::string const &info) {
  TimingScope<3> timing_scope("baz", "");
}

void bar(std::string const &info) {
  TimingScope<2> timing_scope("bar", info);

  baz("");
}

void foo(std::string const &info) {
  TimingScope<2> timing_scope("foo", info);

  bar(info);
  baz("");
}

int main(std::string const &info) {
  TimingScope<1> timing_scope("main", "");

  foo("C30");
  foo("C20V");
  bar("C1");
  bar("");
}
