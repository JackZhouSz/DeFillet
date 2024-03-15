//
// Created by 13900K on 2024/3/15.
//

#include <CLI/CLI.hpp>


int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string param_file;
        std::string output_dir;
    } args;

    CLI::App app{"DEFILLET Command Line"};
    std::string input_path;
    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-p,--parameter", args.param_file, "Parameter")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    CLI11_PARSE(app, argc, argv);
    return 0;
}