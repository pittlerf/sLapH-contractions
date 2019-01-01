#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

struct QuarklineSpec {
  std::string name;
  ssize_t q1;
  ssize_t q2;

  bool is_loop() const { return q1 == q2; }
};

using TraceSpec = std::vector<QuarklineSpec>;
using TraceSpecs = std::vector<TraceSpec>;

// The following data structure holds the source and sink vertex indices in
// the diagrams. The first part of the pair are the source vertices, the
// second part are the sink vertices.
using Vertices = std::pair<std::vector<int>, std::vector<int>>;

struct DiagramSpec {
  Vertices vertices;
  TraceSpecs traces;
};

using DiagramSpecs = std::map<std::string, DiagramSpec>;

extern DiagramSpecs const diagram_specs;
