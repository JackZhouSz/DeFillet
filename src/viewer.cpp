/********************************************************************
 * Copyright (C) 2015 Liangliang Nan <liangliang.nan@gmail.com>
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++ library
 *      for processing and rendering 3D data.
 *      Journal of Open Source Software, 6(64), 3255, 2021.
 * ------------------------------------------------------------------
 *
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 ********************************************************************/

#include "viewer.h"
#include "easy3d/util/dialogs.h"

#include <iostream>

#include <easy3d/util/file_system.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/renderer/text_renderer.h>
#include <easy3d/renderer/camera.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/manipulator.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/renderer/opengl_util.h>
#include <easy3d/renderer/opengl_error.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/drawable_triangles.h>


#include <3rd_party/imgui/misc/fonts/imgui_fonts_droid_sans.h>
#include <3rd_party/imgui/imgui.h>
#include <3rd_party/imgui/backends/imgui_impl_glfw.h>
#include <3rd_party/imgui/backends/imgui_impl_opengl3.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>


namespace easy3d {

    ImGuiContext* ViewerImGui::context_ = nullptr;

    ViewerImGui::ViewerImGui(
            const std::string& title /* = "Easy3D ImGui Viewer" */,
            int samples /* = 4 */,
            int gl_major /* = 3 */,
            int gl_minor /* = 2 */,
            bool full_screen /* = false */,
            bool resizable /* = true */,
            int depth_bits /* = 24 */,
            int stencil_bits /* = 8 */
    )
            : Viewer(title, samples, gl_major, gl_minor, full_screen, resizable, depth_bits, stencil_bits)
            , alpha_(0.8f)
            , movable_(true)
            , mesh(nullptr)
            , sites(nullptr)
            , vertices(nullptr)
    {
        camera()->setUpVector(vec3(0, 1, 0));
        camera()->setViewDirection(vec3(0, 0, -1));
        camera_->showEntireScene();
    }


    void ViewerImGui::init() {
        Viewer::init();

        if (!context_) {
            // Setup ImGui binding
            IMGUI_CHECKVERSION();

            context_ = ImGui::CreateContext();

            const char* glsl_version = "#version 150";
            ImGui_ImplGlfw_InitForOpenGL(window_, true);
            ImGui_ImplOpenGL3_Init(glsl_version);
            ImGuiIO& io = ImGui::GetIO();
            io.WantCaptureKeyboard = true;
            io.WantTextInput = true;
            io.IniFilename = nullptr;
            ImGui::StyleColorsDark();
            ImGuiStyle& style = ImGui::GetStyle();
            style.FrameRounding = 5.0f;

            // load font
            reload_font();
        }
    }


    double ViewerImGui::pixel_ratio() {
        // Computes pixel ratio for hidpi devices
        int fbo_size[2], win_size[2];
        glfwGetFramebufferSize(window_, &fbo_size[0], &fbo_size[1]);
        glfwGetWindowSize(window_, &win_size[0], &win_size[1]);
        return static_cast<double>(fbo_size[0]) / static_cast<double>(win_size[0]);
    }


    void ViewerImGui::reload_font(int font_size)
    {
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts->Clear();
        io.Fonts->AddFontFromMemoryCompressedTTF(droid_sans_compressed_data, droid_sans_compressed_size, font_size * dpi_scaling());
        io.FontGlobalScale = 1.0f / pixel_ratio();
        ImGui_ImplOpenGL3_DestroyDeviceObjects();
    }


    void ViewerImGui::post_resize(int w, int h) {
        Viewer::post_resize(w, h);
        if (context_) {
            ImGui::GetIO().DisplaySize.x = float(w);
            ImGui::GetIO().DisplaySize.y = float(h);
        }
    }


    bool ViewerImGui::callback_event_cursor_pos(double x, double y) {
        if (ImGui::GetIO().WantCaptureMouse)
            return true;
        else
            return Viewer::callback_event_cursor_pos(x, y);
    }


    bool ViewerImGui::callback_event_mouse_button(int button, int action, int modifiers) {
        if (ImGui::GetIO().WantCaptureMouse)
            return true;
        else
            return Viewer::callback_event_mouse_button(button, action, modifiers);
    }


    bool ViewerImGui::callback_event_keyboard(int key, int action, int modifiers) {
        if (ImGui::GetIO().WantCaptureKeyboard)
            return true;
        else
            return Viewer::callback_event_keyboard(key, action, modifiers);
    }


    bool ViewerImGui::callback_event_character(unsigned int codepoint) {
        if (ImGui::GetIO().WantCaptureKeyboard)
            return true;
        else
            return Viewer::callback_event_character(codepoint);
    }


    bool ViewerImGui::callback_event_scroll(double dx, double dy) {
        if (ImGui::GetIO().WantCaptureMouse)
            return true;
        else
            return Viewer::callback_event_scroll(dx, dy);
    }


    void ViewerImGui::cleanup() {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();

        ImGui::DestroyContext(context_);

        Viewer::cleanup();
    }


    void ViewerImGui::pre_draw() {
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        Viewer::pre_draw();
    }

    void ViewerImGui::draw() const {
        if(mesh && mesh->renderer()->is_visible()) {
            std::size_t count = 0;
            for (auto d : mesh->renderer()->lines_drawables()) {
                if (d->is_visible()) {
                    d->draw(camera()); easy3d_debug_log_gl_error;
                    ++count;
                }
            }

            for (auto d : mesh->renderer()->points_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }

            if (count > 0) {
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(0.5f, -0.0001f);
            }
            for (auto d : mesh->renderer()->triangles_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }
            if (count > 0)
                glDisable(GL_POLYGON_OFFSET_FILL);
        }

        if(sites && sites->renderer()->is_visible()) {
            std::size_t count = 0;
            for (auto d : sites->renderer()->lines_drawables()) {
                if (d->is_visible()) {
                    d->draw(camera()); easy3d_debug_log_gl_error;
                    ++count;
                }
            }

            for (auto d : sites->renderer()->points_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }

            if (count > 0) {
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(0.5f, -0.0001f);
            }
            for (auto d : sites->renderer()->triangles_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }
            if (count > 0)
                glDisable(GL_POLYGON_OFFSET_FILL);
        }

        if(vertices && vertices->renderer()->is_visible()) {
            std::size_t count = 0;
            for (auto d : vertices->renderer()->lines_drawables()) {
                if (d->is_visible()) {
                    d->draw(camera()); easy3d_debug_log_gl_error;
                    ++count;
                }
            }

            for (auto d : vertices->renderer()->points_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }

            if (count > 0) {
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(0.5f, -0.0001f);
            }
            for (auto d : vertices->renderer()->triangles_drawables()) {
                if (d->is_visible())
                    d->draw(camera()); easy3d_debug_log_gl_error;
            }
            if (count > 0)
                glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }
    void ViewerImGui::post_draw() {
        draw_dashboard();
        draw_corner_axes();
    }



    void ViewerImGui::draw_dashboard(){
        static bool my_tool_active = true;
        int w, h;
        glfwGetWindowSize(window_, &w, &h);
        ImGui::SetNextWindowSize(ImVec2(400, h), ImGuiCond_Always);
        ImGui::SetNextWindowPos(ImVec2(w - 400, 0), ImGuiCond_Always);
        ImGui::Begin("Dashboard", &my_tool_active, ImGuiWindowFlags_MenuBar);
        ImGui::SetNextWindowContentSize(ImVec2(0, 0));
        if (ImGui::CollapsingHeader("Open/Save", ImGuiTreeNodeFlags_DefaultOpen)) {

            if (ImGui::Button("Open", ImVec2(200, 30))) {
                const std::string title("Please choose a file");
                const std::string &default_path = "../data/";
                const std::vector<std::string> &filters = {
                        "Surface Mesh (*.obj *.ply)", "*.obj *.ply"
                };
                const std::vector<std::string> &file_names = easy3d::dialog::open(title, default_path, filters, true);
                if(mesh) {
                    delete mesh->renderer();
                    delete mesh->manipulator();
                    delete mesh;
                }
                mesh = easy3d::SurfaceMeshIO::load(file_names[0]);
                auto renderer = new easy3d::Renderer(mesh, true);
                mesh->set_renderer(renderer);
                auto manipulator = new Manipulator(mesh);
                mesh->set_manipulator(manipulator);
                Viewer::fit_screen(mesh);

            }
            ImGui::SameLine();
            if (ImGui::Button("Save", ImVec2(200, 30))) {

            }
        }
        if (ImGui::CollapsingHeader("Visible", ImGuiTreeNodeFlags_DefaultOpen)) {
            static bool show_mesh = true;
            ImGui::Checkbox("Panel Movable", &show_mesh);
            if(mesh && mesh->renderer()) {
                mesh->renderer()->set_visible(show_mesh);
            }
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    }
}
