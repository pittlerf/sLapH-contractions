#include <timings.hpp>

void foo();
void bar();
void baz();

void foo() {
  TimingScope<2> timing_scope("foo", "");

  bar();
  baz();
}

void bar() {
  TimingScope<2> timing_scope("bar", "");

  baz();
}

void baz() {
  TimingScope<3> timing_scope("baz", "");
}

int main() {
  TimingScope<1> timing_scope("main", "");

  foo();
  bar();
}
