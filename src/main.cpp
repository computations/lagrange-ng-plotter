/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2022 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include <genesis/genesis.hpp>
#include <numeric>
#include <string>

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <limits>
#include <nlohmann/json.hpp>
#include <unordered_map>

using namespace std;

struct CLIOptions {
  string results_json_filename;
  string output_filename;
};

struct StateResult {
  std::string distriution_name;
  size_t      distribution;
  double      ratio;
};

struct PieChartEntry {
  double value;
  double color_index;
  string entry_label;
};

struct NodeResult {
  std::string           node_id;
  vector<StateResult>   state_results;
  vector<PieChartEntry> pie_chart_entries;
};

struct Attributes {
  size_t regions;
  size_t state_count;
  size_t max_areas;
  string newick_tree;
};

// Stolen code from lagrange-ng to compute the next dist
inline auto next_dist(size_t d, uint32_t n) -> size_t {
  d += 1;
  while (static_cast<size_t>(__builtin_popcountll(d)) > n) { d++; }
  return d;
}

std::vector<double>
get_node_values(const unordered_map<string, NodeResult> &node_results,
                const string                            &node_id) {
  const auto &nr = node_results.at(node_id);

  vector<double> values;

  for (const auto &pce : nr.pie_chart_entries) { values.push_back(pce.value); }

  return values;
}

SvgGroup make_random_pie_chart(std::vector<Color> const &colors) {
  // Super simple, no need for C++ random classes for now.
  std::vector<double> values;
  size_t const        n = 2 + rand() % 8;
  for (size_t i = 0; i < n; ++i) { values.emplace_back(rand() % 10); }

  // Make the pie chart.
  return make_svg_pie_chart(values, colors, 10.0);
}

SvgDocument
make_pie_chart_tree(Tree                                    &tree,
                    unordered_map<string, NodeResult> const &node_results,
                    vector<string> const                    &node_map,
                    LayoutParameters const                  &params,
                    std::vector<Color> const                &colors) {
  // Make a layout tree. We need a pointer to it in order to allow for the two
  // different classes (circular/rectancular) to be returned here. Make it a
  // unique ptr for automatic cleanup.
  auto layout = [&]() -> std::unique_ptr<LayoutBase> {
    if (params.shape == LayoutShape::kCircular) {
      return utils::make_unique<CircularLayout>(
          tree, params.type, params.ladderize);
    }
    if (params.shape == LayoutShape::kRectangular) {
      return utils::make_unique<RectangularLayout>(
          tree, params.type, params.ladderize);
    }
    throw std::runtime_error("Unknown Tree shape parameter.");
  }();

  // Add pie chars at the leaf nodes.
  auto node_shapes = std::vector<utils::SvgGroup>(layout->tree().node_count());
  for (size_t i = 0; i < layout->tree().node_count(); ++i) {
    if (is_leaf(layout->tree().node_at(i))) { continue; }

    node_shapes[i] = make_svg_pie_chart(
        get_node_values(node_results, node_map[i]), colors, 10.0);
  }
  layout->set_node_shapes(node_shapes);
  layout->set_edge_strokes(params.stroke);

  for (auto &node : tree.nodes()) {
    if (!is_inner(node)) { continue; }
    node.data<CommonNodeData>().name = "";
  }

  return layout->to_svg_document();
}

void add_color_list_legend(SvgDocument               &svg_doc,
                           LayoutParameters const    &params,
                           std::vector<Color> const  &colors,
                           std::vector<string> const &color_labels) {
  // Make the color list. There also is a version for gradients
  // make_svg_color_bar(), but we want a discrete list here instead for the pie
  // chart.
  auto svg_color_list = make_svg_color_list(colors, color_labels);

  // Height, as fraction of overall tree height.
  double const legend_height = 0.25;

  // Move it to the bottom right corner.
  if (params.shape == LayoutShape::kCircular) {
    // Circular trees are centered around 0, so we want to start from there.
    svg_color_list.transform.append(utils::SvgTransform::Translate(
        1.5 * svg_doc.bounding_box().width() / 2.0,
        (0.5 - legend_height) * svg_doc.bounding_box().height()));
  }
  if (params.shape == LayoutShape::kRectangular) {
    // Rect trees don't have a center, so we need a different approach,
    // and just treat the whole tree as a box, where 0,0 is the bottom left.
    svg_color_list.transform.append(utils::SvgTransform::Translate(
        1.1 * svg_doc.bounding_box().width(),
        (1.0 - legend_height) * svg_doc.bounding_box().height()));
  }

  // Apply a scale factor that scales the box to be a fraction of the figure
  // height. The denominator is the number items in the list times their height
  // (hardcoded 15px, used by make_svg_color_list - sorry for that...)
  auto const sf = ((legend_height * svg_doc.bounding_box().height())
                   / (static_cast<double>(colors.size()) * 15.0));
  svg_color_list.transform.append(utils::SvgTransform::Scale(sf));

  // Add it to the svg doc.
  svg_doc.add(svg_color_list);
  svg_doc.margin.left  = 10;
  svg_doc.margin.right = 50;
}

CLIOptions parse_options(int argc, char **argv) {
  CLIOptions           options;
  static struct option long_options[] = {{"results", required_argument, 0, 0},
                                         {"output", required_argument, 0, 0},
                                         {0, 0, 0, 0}};
  int                  c              = 0;
  int                  option_index   = 0;
  optind                              = 0;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index))
         == 0) {
    switch (option_index) {
    case 0: // results
      options.results_json_filename = optarg;
      break;
    case 1: // output
      options.output_filename = optarg;
      break;
    case '?':
    case ':':
      std::exit(1);
      break;
    default:
      throw std::runtime_error{"CLI Argument was not recognized"};
    }
  }
  return options;
}

void check_options(const CLIOptions &options) {
  if (options.results_json_filename.empty()) {
    LOG_INFO << "A results file is required, please provide a results file "
                "with the --results switch";
    std::exit(1);
  }

  if (options.output_filename.empty()) {
    LOG_INFO << "An output filename is required, please provide an output "
                "filename with the --output switch";
    std::exit(1);
  }
}

NodeResult parse_node_result(const nlohmann::json &json) {
  NodeResult node_result;
  node_result.node_id = json["number"];
  for (const auto &el : json["states"]) {
    StateResult sr;
    sr.distribution     = el["distribution"];
    sr.distriution_name = [](vector<string> ranges) -> std::string {
      if (ranges.empty()) { return {}; }
      if (ranges.size() == 1) { return ranges.front(); }
      auto fold_operation
          = [](std::string a, std::string b) { return a + ", " + b; };
      return std::accumulate(
          ranges.begin() + 1, ranges.end(), ranges.front(), fold_operation);
    }(el["regions"].template get<vector<string>>());
    sr.ratio = el["ratio"];
    node_result.state_results.push_back(sr);
  }
  return node_result;
}

unordered_map<string, NodeResult>
get_node_results_from_json(const nlohmann::json &json) {
  unordered_map<string, NodeResult> node_results;

  for (const auto &el : json["node-results"]) {
    node_results[el["number"]] = parse_node_result(el);
  }

  return node_results;
}

void select_pie_chart_results(unordered_map<string, NodeResult> &node_results,
                              size_t                             states,
                              size_t                             max_areas,
                              size_t                             color_count) {
  std::vector<std::pair<size_t, double>> summed_values;

  for (size_t i = 0, dist = 0; i < states;
       i++, dist          = next_dist(dist, max_areas)) {
    summed_values.emplace_back(dist, 0.0);
  }

  for (const auto &n : node_results) {
    auto node_results = n.second;
    for (size_t index = 0; index < node_results.state_results.size(); ++index) {
      const auto &sr = node_results.state_results[index];
      assert(summed_values[index].first == sr.distribution);
      summed_values[index].second += sr.ratio;
    }
  }

  std::sort(summed_values.begin(),
            summed_values.end(),
            [](auto a, auto b) -> bool { return a.second > b.second; });

  summed_values.resize(color_count - 1);

  for (auto &n : node_results) {
    auto &node_result = n.second;

    double total = 1.0;
    for (size_t j = 0; j < summed_values.size(); ++j) {
      auto          sr = node_result.state_results[summed_values[j].first];
      PieChartEntry pce;
      pce.color_index  = j;
      pce.value        = sr.ratio;
      pce.entry_label  = sr.distriution_name;
      total           -= sr.ratio;
      node_result.pie_chart_entries.push_back(pce);
    }

    assert(total >= -std::numeric_limits<double>::epsilon());
    if (total < 0.0) { total = 0.0; }

    PieChartEntry pce_other;
    pce_other.value       = total;
    pce_other.color_index = color_count - 1;
    pce_other.entry_label = "Other";
    node_result.pie_chart_entries.push_back(pce_other);
  }
}

Attributes get_attributes(const nlohmann::json &json) {
  Attributes at;

  at.newick_tree = json["attributes"]["nodes-tree"];
  at.regions     = json["attributes"]["regions"];
  at.state_count = json["attributes"]["state-count"];
  at.max_areas   = json["attributes"]["max-areas"];

  return at;
}

vector<string> produce_node_map(const Tree &tree) {
  vector<string> node_map;
  node_map.resize(tree.node_count());
  for (size_t node_index = 0; node_index < tree.node_count(); node_index++) {
    const auto &node = tree.node_at(node_index);
    if (is_leaf(node)) { continue; }
    node_map[node_index] = node.data<CommonNodeData>().name;
  }
  return node_map;
}

int main(int argc, char **argv) {
  // Activate logging.
  utils::Logging::log_to_stdout();
  utils::Logging::details.time = true;
  LOG_INFO << "Started";

  // For testing, we use one of our default color schemes for categorical data.
  auto const colors = color_list_set2();

  auto options = parse_options(argc, argv);

  fstream        json_file(options.results_json_filename);
  nlohmann::json results_json = nlohmann::json::parse(json_file);

  auto attributes = get_attributes(results_json);
  LOG_INFO << "Parsed attributes";
  auto node_results = get_node_results_from_json(results_json);
  LOG_INFO << "Parsed results";

  select_pie_chart_results(node_results,
                           attributes.state_count,
                           attributes.max_areas,
                           colors.size());
  LOG_INFO << "Finished precomputation";

  // Read the tree.
  auto tree
      = CommonTreeNewickReader().read(from_string(attributes.newick_tree));
  LOG_INFO << "Found tree with " << leaf_node_count(tree)
           << " tips in tree file";

  auto node_map = produce_node_map(tree);
  LOG_INFO << "Finished node map";

  // clear the names for inner nodes of the tree
  for (auto &node : tree.nodes()) {
    if (is_leaf(node)) { continue; }
    node.data<CommonNodeData>().name.clear();
  }

  // Set the tree params as needed.
  LayoutParameters params;
  // params.shape = LayoutShape::kCircular;
  params.shape           = LayoutShape::kRectangular;
  params.type            = LayoutType::kCladogram;
  params.ladderize       = true;
  params.stroke          = SvgStroke(Color(), 1.0);
  params.stroke.line_cap = utils::SvgStroke::LineCap::kRound;

  // Get the layout as an svg doc.
  auto svg_doc
      = make_pie_chart_tree(tree, node_results, node_map, params, colors);

  LOG_INFO << "Made SVG";
  vector<string> color_labels;

  for (const auto &pce : node_results.begin()->second.pie_chart_entries) {
    color_labels.push_back(pce.entry_label);
  }
  LOG_INFO << "Made labels";

  // Add a legend for our colors.
  add_color_list_legend(svg_doc, params, colors, color_labels);

  // Write the whole svg doc to file.
  std::ofstream ofs;
  utils::file_output_stream(options.output_filename, ofs);
  svg_doc.write(ofs);

  LOG_INFO << "Finished";
  return 0;
}
