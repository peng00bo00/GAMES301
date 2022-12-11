#include "tutte.h"
#include "parametrization.h"

#include <algorithm>
#include <igl/readOBJ.h>
#include <igl/Timer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <matplot/matplot.h>

using rt = igl::opengl::ViewerCore::RotationType;

int main()
{
    Eigen::MatrixXd V, V_uv;
    Eigen::MatrixXi F;

    ParamData data;

    igl::opengl::glfw::Viewer viewer;
    igl::Timer timer;

    auto figure = matplot::figure(true);
    figure->size(1280, 960);

    unsigned int left_view, right_view;
    const auto mesh_id = viewer.data_list[0].id;
    const auto uv_id = viewer.append_mesh(false);

    int current_method = static_cast<int>(weight_method::harmonic);
    int current_border = static_cast<int>(border_type::circle);
    int current_energy = static_cast<int>(energy_type::sd);
    int current_max_iter = 3000;
    int animation_frame_count = 999999;
    bool animation = false;
    std::vector<double> animation_figure;

    double uv_scale_param;

    viewer.callback_init = [&](igl::opengl::glfw::Viewer&)
    {
        viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
        left_view = viewer.core_list[0].id;
        right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));

        viewer.data(mesh_id).set_visible(false, right_view);

        viewer.core(left_view).set_rotation_type(rt::ROTATION_TYPE_TRACKBALL);
        viewer.core(right_view).set_rotation_type(rt::ROTATION_TYPE_NO_ROTATION);

        viewer.core(right_view).orthographic = true;
        viewer.core(right_view).background_color = {255, 255, 255, 255};

        return false;
    };
    viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer& v, int w, int h)
    {
        v.core(left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
        v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
        return true;
    };

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Load##Mesh", ImVec2(-1, 0)))
            {
                std::string fname = igl::file_dialog_open();

                if (fname.length() != 0)
                {
                    viewer.data(mesh_id).clear();
                    viewer.data(uv_id).clear();

                    igl::read_triangle_mesh(fname.c_str(), V, F);

                    tutte(V, F, V_uv, current_method, current_border);

                    data = projected_newton_precompute(V, F);

                    uv_scale_param = 15 * (1. / std::sqrt(data.mesh_area));

                    viewer.data(mesh_id).set_mesh(V, F);
                    viewer.data(mesh_id).set_uv(V_uv * 5);

                    viewer.data(uv_id).set_mesh(V_uv, F);

                    viewer.core(left_view).align_camera_center(V, F);
                    viewer.core(right_view).align_camera_center(V_uv, F);

                    viewer.data(mesh_id).show_lines = false;
                    viewer.data(mesh_id).show_texture = true;
                }
            }
        }

        if (ImGui::CollapsingHeader("Viewer Settings", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::DragFloat("UV Line Width", &viewer.data(uv_id).line_width, 0.01f, 0.1f, 5.0f, "%.2f",
                             ImGuiSliderFlags_AlwaysClamp);
        }

        if (ImGui::CollapsingHeader("Parametrization Settings", ImGuiTreeNodeFlags_DefaultOpen))
        {
            const char* energy_items[] = {"Symmetric Dirichlet", "ARAP", "MIPS"};

            ImGui::Text("Distortion Energy Type: ");

            ImGui::PushItemWidth(1.2f * ImGui::CalcTextSize(energy_items[0]).x + ImGui::GetStyle().FramePadding.x);
            ImGui::Combo("", &current_energy, energy_items, IM_ARRAYSIZE(energy_items));
            ImGui::PopItemWidth();

            ImGui::Text("Max iteration number");
            ImGui::SameLine(0, ImGui::GetStyle().FramePadding.x * 2.0f);
            ImGui::PushItemWidth((ImGui::GetContentRegionAvail().x - ImGui::GetStyle().FramePadding.x) * 0.75f);
            ImGui::DragInt("", &current_max_iter, 1, 1, 10000, "%d", ImGuiSliderFlags_AlwaysClamp);
            ImGui::PopItemWidth();

            if (ImGui::Checkbox("Animation", &animation))
            {
                if (animation)
                {
                    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& v)
                    {
                        if (animation && animation_frame_count < current_max_iter)
                        {
                            if(animation_frame_count == 0)
                            {
                                tutte(V, F, V_uv, current_method, current_border);
                                v.data(mesh_id).set_uv(V_uv * 5);
                                v.data(uv_id).clear();
                                v.data(uv_id).set_mesh(V_uv, F);
                            }

                            const auto& [old_energy, new_energy] = projected_newton_solver_animation(
                                V_uv, V, F, data, current_energy);
                            if (animation_frame_count == 0)
                            {
                                animation_figure.push_back(old_energy);
                            }

                            if (new_energy == -1)
                            {
                                animation_frame_count = 99999;
                            }
                            else
                            {
                                animation_figure.push_back(new_energy);
                            }

                            v.data(mesh_id).set_uv(V_uv * uv_scale_param);
                            v.data(uv_id).clear();
                            v.data(uv_id).set_mesh(V_uv, F);

                            viewer.core(left_view).align_camera_center(V, F);
                            viewer.core(right_view).align_camera_center(V_uv, F);

                            animation_frame_count++;

                            figure->current_axes()->clear();
                            figure->current_axes()->plot(animation_figure)->line_width(2);
                            figure->current_axes()->xlabel("Iteration");
                            figure->current_axes()->ylabel("Energy");
                            figure->draw();
                        }
                        return false;
                    };
                }
                else
                {
                    animation_frame_count = 99999;
                    animation_figure = {};

                    viewer.callback_pre_draw = nullptr;
                    viewer.core(left_view).is_animating = false;
                    viewer.core(right_view).is_animating = false;
                }
            }
        }

        if (ImGui::Button("Run", ImVec2(-1, 0)))
        {
            if (tutte(V, F, V_uv, current_method, current_border))
            {
                if (animation)
                {
                    animation_frame_count = 0;
                    animation_figure = {};

                    viewer.core(left_view).is_animating = true;
                    viewer.core(right_view).is_animating = true;
                }
                else
                {
                    std::cout << "Start parametrization.\n";
                    timer.start();

                    const auto& energy_vector = projected_newton_solver(V_uv, V, F, data, current_energy,
                                                                        current_max_iter);

                    timer.stop();
                    std::cout << "Done.\nTakes " << timer.getElapsedTime() << " seconds.\n";

                    const auto last_idx = std::find(energy_vector.begin(), energy_vector.end(), -1);

                    std::cout << "Init Energy: " << energy_vector[0] << " Final Energy: " << (
                            last_idx != energy_vector.end() ? *(last_idx - 1) : energy_vector.back()) << ". It takes "
                        <<
                        std::distance(energy_vector.begin(), last_idx) - 1 << " iteration to finish." << std::endl;

                    viewer.data(mesh_id).set_uv(V_uv * uv_scale_param);

                    viewer.data(uv_id).set_vertices(V_uv);

                    viewer.core(left_view).align_camera_center(V, F);
                    viewer.core(right_view).align_camera_center(V_uv, F);

                    const std::vector draw_energy(energy_vector.begin(), last_idx);

                    figure->current_axes()->clear();
                    figure->current_axes()->plot(draw_energy)->line_width(2);
                    figure->current_axes()->xlabel("Iteration");
                    figure->current_axes()->ylabel("Energy");
                    figure->draw();
                }
            }
        }
    };

    viewer.data(mesh_id).show_lines = false;
    viewer.data(mesh_id).show_texture = true;

    viewer.data(uv_id).show_lines = true;
    viewer.data(uv_id).show_texture = false;
    viewer.data(uv_id).show_faces = false;

    viewer.core(left_view).is_animating = false;
    viewer.core(right_view).is_animating = false;

    viewer.launch();
}
