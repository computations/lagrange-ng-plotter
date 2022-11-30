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

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <fstream>
#include <getopt.h>

using namespace std;

struct CLIOptions {
  string results_json_filename;
  string output_filename;
};

SvgGroup make_random_pie_chart(std::vector<Color> const &colors) {
  // Super simple, no need for C++ random classes for now.
  std::vector<double> values;
  size_t const        n = 2 + rand() % 8;
  for (size_t i = 0; i < n; ++i) { values.emplace_back(rand() % 10); }

  // Make the pie chart.
  return make_svg_pie_chart(values, colors, 10.0);
}

SvgDocument make_pie_chart_tree(Tree const               &tree,
                                LayoutParameters const   &params,
                                std::vector<Color> const &colors) {
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
    node_shapes[i] = make_random_pie_chart(colors);
  }
  layout->set_node_shapes(node_shapes);
  // layout->set_edge_strokes(params.stroke);

  return layout->to_svg_document();
}

void add_color_list_legend(SvgDocument              &svg_doc,
                           LayoutParameters const   &params,
                           std::vector<Color> const &colors) {
  // For testing, we just make a fake list of labels, numerbing the colors...
  std::vector<std::string> color_labels;
  for (size_t i = 0; i < colors.size(); ++i) {
    color_labels.push_back("Color " + std::to_string(i + 1));
  }

  // Make the color list. There also is a version for gradients
  // make_svg_color_bar(), but we want a discrete list here instead for the pie
  // chart.
  auto svg_color_list = make_svg_color_list(colors, color_labels);

  // Height, as fraction of overall tree height.
  double const legend_height = 0.5;

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
        1.5 * svg_doc.bounding_box().width(),
        (1.0 - legend_height) * svg_doc.bounding_box().height()));
  }

  // Apply a scale factor that scales the box to be a fraction of the figure
  // height. The denominator is the number items in the list times their height
  // (hardcoded 15px, used by make_svg_color_list - sorry for that...)
  auto const sf = ((legend_height * svg_doc.bounding_box().height()) /
                   (static_cast<double>(colors.size()) * 15.0));
  svg_color_list.transform.append(utils::SvgTransform::Scale(sf));

  // Add it to the svg doc.
  svg_doc.add(svg_color_list);
  svg_doc.margin.right += 100;
}

CLIOptions parse_options(int argc, char **argv) {
  CLIOptions           options;
  static struct option long_options[] = {{"results", required_argument, 0, 0},
                                         {"output", required_argument, 0, 0},
                                         {0, 0, 0, 0}};
  int                  c              = 0;
  int                  option_index   = 0;
  optind                              = 0;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) ==
         0) {
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

nlohmann::json parse_results_json(const string &filename) {
  fstream json_file(filename);
  return nlohmann::json::parse(json_file);
}

string get_tree_newick(const nlohmann::json &results_json) {
  return results_json["attributes"]["nodes-tree"];
}

int main(int argc, char **argv) {
  // Activate logging.
  utils::Logging::log_to_stdout();
  utils::Logging::details.time = true;
  LOG_INFO << "Started";

  auto options      = parse_options(argc, argv);
  auto results_json = parse_results_json(options.results_json_filename);

  // Read the tree.
  auto const tree =
      CommonTreeNewickReader().read(from_string(get_tree_newick(results_json)));
  LOG_INFO << "Found tree with " << leaf_node_count(tree)
           << " tips in tree file";

  // Set the tree params as needed.
  LayoutParameters params;
  // params.shape = LayoutShape::kCircular;
  params.shape           = LayoutShape::kRectangular;
  params.type            = LayoutType::kCladogram;
  params.ladderize       = true;
  params.stroke          = SvgStroke(Color(), 1.0);
  params.stroke.line_cap = utils::SvgStroke::LineCap::kRound;

  // For testing, we use one of our default color schemes for categorical data.
  auto const colors = color_list_set1();

  // Get the layout as an svg doc.
  auto svg_doc = make_pie_chart_tree(tree, params, colors);

  // Add a legend for our colors.
  add_color_list_legend(svg_doc, params, colors);

  // Write the whole svg doc to file.
  std::ofstream ofs;
  utils::file_output_stream(options.output_filename, ofs);
  svg_doc.write(ofs);

  LOG_INFO << "Finished";
  return 0;
}
