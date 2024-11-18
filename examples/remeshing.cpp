// Copyright 2011-2019 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include <pmp/visualization/mesh_viewer.h>
#include <pmp/algorithms/features.h>
#include <pmp/algorithms/remeshing.h>
#include <pmp/algorithms/utilities.h>

#include <imgui.h>

#include "pmp/algorithms/utilities.h"

#include "pmp/io/io.h"

using namespace pmp;

class Viewer : public MeshViewer
{
public:
    Viewer(const char* title, int width, int height);

protected:
    void process_imgui() override;
};

void mesh_statistics(SurfaceMesh& mesh)
{
    float avg_edge_len = 0;
    float max_edge_len = std::numeric_limits<float>::min();
    float min_edge_len = std::numeric_limits<float>::max();
    float avg_vertex_valence = 0;
    float max_vertex_valence = std::numeric_limits<int>::min();
    float min_vertex_valence = std::numeric_limits<int>::max();

    for (auto e : mesh.edges())
    {
        auto v0 = mesh.vertex(e, 0);
        auto v1 = mesh.vertex(e, 1);

        float l = distance(mesh.vertex_property<Point>("v:point")[v0],
                           mesh.vertex_property<Point>("v:point")[v1]);

        avg_edge_len += l;
        max_edge_len = std::max(max_edge_len, l);
        min_edge_len = std::min(min_edge_len, l);
    }

    avg_edge_len /= mesh.n_edges();

    for (auto v : mesh.vertices())
    {
        float valence = float(mesh.valence(v));

        avg_vertex_valence += valence;
        max_vertex_valence = std::max(max_vertex_valence, valence);
        min_vertex_valence = std::min(min_vertex_valence, valence);
    }
    avg_vertex_valence /= mesh.n_vertices();

    printf(
        "\n Mesh Stats: Avg Edge Length= %f, Max Edge Length= %f, Min Edge "
        "Length= %f, Avg Vertex Valence= %f, Max Vertex Valence= %f, Min "
        "Vertex Valence= %f\n",
        avg_edge_len, max_edge_len, min_edge_len, avg_vertex_valence,
        max_vertex_valence, min_vertex_valence);
}

Viewer::Viewer(const char* title, int width, int height)
    : MeshViewer(title, width, height)
{
    set_draw_mode("Hidden Line");
    crease_angle_ = 0.0;
}

void Viewer::process_imgui()
{
    MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Features", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();

        static int feature_angle = 70;
        ImGui::PushItemWidth(80);
        ImGui::SliderInt("##feature_angle", &feature_angle, 1, 180);
        ImGui::PopItemWidth();
        ImGui::SameLine();
        if (ImGui::Button("Detect Features"))
        {
            clear_features(mesh_);
            detect_features(mesh_, feature_angle);
            update_mesh();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Uniform Remeshing",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();

        ImGui::PushItemWidth(80);

        static double edge_length{0.01};
        ImGui::InputDouble("Edge Length", &edge_length, 0, 0, "%g");
        ImGui::SameLine();
        if (ImGui::Button("Mean"))
        {
            edge_length = mean_edge_length(mesh_);
        }

        static int n_iterations{10};
        ImGui::SliderInt("Iterations##uniform", &n_iterations, 1, 20);

        static bool use_projection{false};
        ImGui::Checkbox("Use Projection##uniform", &use_projection);

        static bool scale_lengths{false};
        ImGui::Checkbox("Scale Lengths##uniform", &scale_lengths);

        ImGui::Spacing();

        if (ImGui::Button("Remesh##uniform"))
        {
            try
            {
                auto scaling = scale_lengths ? bounds(mesh_).size() : 1.0;
                uniform_remeshing(mesh_, edge_length * scaling, n_iterations,
                                  use_projection);
            }
            catch (const InvalidInputException& e)
            {
                std::cerr << e.what() << std::endl;
                return;
            }
            update_mesh();
        }

        ImGui::PopItemWidth();
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Adaptive Remeshing",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();

        ImGui::PushItemWidth(80);

        static double min_length{0.001};
        ImGui::InputDouble("Min. Edge Length", &min_length, 0, 0, "%g");

        static double max_length{0.05};
        ImGui::InputDouble("Max. Edge Length", &max_length, 0, 0, "%g");

        static double max_error{0.0005};
        ImGui::InputDouble("Max. Error", &max_error, 0, 0, "%g");

        static int n_iterations{10};
        ImGui::SliderInt("Iterations##adaptive", &n_iterations, 1, 20);

        static bool use_projection{true};
        ImGui::Checkbox("Use Projection##adaptive", &use_projection);

        static bool scale_lengths{true};
        ImGui::Checkbox("Scale Lengths##adaptive", &scale_lengths);

        ImGui::Spacing();

        if (ImGui::Button("Remesh##adaptive"))
        {
            auto scaling = scale_lengths ? bounds(mesh_).size() : 1.0;
            try
            {
                adaptive_remeshing(mesh_, min_length * scaling,
                                   max_length * scaling, max_error * scaling,
                                   n_iterations, use_projection);
            }
            catch (const InvalidInputException& e)
            {
                std::cerr << e.what() << std::endl;
                return;
            }
            update_mesh();
        }

        ImGui::PopItemWidth();
    }
}

int main(int argc, char** argv)
{
    SurfaceMesh mesh_;

    read(mesh_, std::string(argv[1]));

    int n_iterations = 1;

    double edge_length = mean_edge_length(mesh_);

    std::cout << "\n #### Input #Faces= " << mesh_.n_faces()
              << " #Vertices=" << mesh_.n_vertices() << "\n";
    mesh_statistics(mesh_);

    auto start = std::chrono::high_resolution_clock::now();

    uniform_remeshing(mesh_, edge_length, n_iterations, false);

    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "\n $$$$ Uniform remeshing took = "
              << std::chrono::duration<float, std::milli>(stop - start).count()
              << "(ms)\n";

    std::cout << "\n #### Output #Faces= " << mesh_.n_faces()
              << " #Vertices=" << mesh_.n_vertices() << "\n";
    mesh_statistics(mesh_);

    write(mesh_, std::string("pmp_remesh.obj"));

    //#ifndef __EMSCRIPTEN__
    //    Viewer window("Remeshing", 800, 600);
    //    if (argc == 2)
    //        window.load_mesh(argv[1]);
    //    return window.run();
    //#else
    //    Viewer window("Remeshing", 800, 600);
    //    window.load_mesh(argc == 2 ? argv[1] : "input.off");
    //    return window.run();
    //#endif
}
