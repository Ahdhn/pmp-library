// Copyright 2011-2019 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include <pmp/visualization/mesh_viewer.h>
#include <pmp/algorithms/sec.h>
#include <imgui.h>

using namespace pmp;

class Viewer : public MeshViewer
{
public:
    Viewer(const char* title, int width, int height);

protected:
    void process_imgui() override;
};

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

    if (ImGui::CollapsingHeader("Decimation", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();
        ImGui::PushItemWidth(80);

        static int target_percentage = 10;
        ImGui::SliderInt("Number of Vertices (%)", &target_percentage, 1, 99);

        if (ImGui::Button("Decimate"))
        {
            try
            {
                auto start = std::chrono::high_resolution_clock::now();
                auto nv = mesh_.n_vertices() * 0.01 * target_percentage;

                shortest_edge_collapse(mesh_, nv);

                auto stop = std::chrono::high_resolution_clock::now();

                std::cout << "\n $$$$ Shortest Edge Collapse took = "
                          << std::chrono::duration<float, std::milli>(stop -
                                                                      start)
                                 .count()
                          << "(ms)\n";
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
#ifndef __EMSCRIPTEN__
    Viewer window("Decimation", 800, 600);
    if (argc == 2)
        window.load_mesh(argv[1]);
    return window.run();
#else
    Viewer window("Decimation", 800, 600);
    window.load_mesh(argc == 2 ? argv[1] : "input.off");
    return window.run();
#endif
}
