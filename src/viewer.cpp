#include "viewer.h"
#include "easy3d/util/dialogs.h"

#include <iostream>
#include <cstdio>
#include <thread>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <filesystem>

#include <easy3d/util/file_system.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/renderer/text_renderer.h>
#include <easy3d/renderer/camera.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/manipulator.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/util/file_system.h>
#include <easy3d/renderer/opengl_util.h>
#include <easy3d/renderer/opengl_error.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/drawable_triangles.h>



#include <easy3d/util/timer.h>

#include <3rd_party/imgui/misc/fonts/imgui_fonts_droid_sans.h>
#include <3rd_party/imgui/imgui.h>
#include <3rd_party/imgui/backends/imgui_impl_glfw.h>
#include <3rd_party/imgui/backends/imgui_impl_opengl3.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>


#include <igl/jet.h>

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
            , fillet_seg(nullptr)
            , show_mesh(false)
            , eps(0.03), s(10), radius(0.06)
            , min_score(0.5), alpha(0.5)
            , std_ratio(0.3), num_sor_iter(3), nb_neighbors(30){
        camera()->setUpVector(vec3(0, 1, 0));
        camera()->setViewDirection(vec3(0, 0, -1));
        camera_->showEntireScene();
        cur_work_dir = easy3d::file_system::current_working_directory();
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

        if(state == UPDATE_SCORING) {
            for(auto m : models_) {
                if(m) {
                    delete_model(m);
                }
            }
            models_.clear();
            mesh = add_model(scoring_mesh_path);
            sites = add_model(scoring_sites_path);
            vertices = add_model(scoring_vertices_path);
            state = NOTHING;
        }



        Viewer::pre_draw();
    }


    void ViewerImGui::post_draw() {
        draw_dashboard();
        draw_corner_axes();
    }


    void ViewerImGui::draw_dashboard(){
        static bool my_tool_active = true;
        int w, h;
        glfwGetWindowSize(window_, &w, &h);
        ImGui::SetNextWindowSize(ImVec2(325, h), ImGuiCond_Always);
        ImGui::SetNextWindowPos(ImVec2(w - 325, 0), ImGuiCond_Always);
        ImGui::Begin("Dashboard", &my_tool_active, ImGuiWindowFlags_MenuBar);
        ImGui::SetNextWindowContentSize(ImVec2(0, 0));
        if (ImGui::CollapsingHeader("Open/Save", ImGuiTreeNodeFlags_DefaultOpen)) {

            if (ImGui::Button("Open", ImVec2(150, 30))) {
                open();
            }
            ImGui::SameLine();
            if (ImGui::Button("Save", ImVec2(150, 30))) {

            }
        }
        ImGui::Separator();
        if (ImGui::CollapsingHeader("Visible", ImGuiTreeNodeFlags_DefaultOpen)) {
            if(mesh && mesh->renderer()) {
                show_mesh = mesh->renderer()->is_visible();
            }
            ImGui::Checkbox("mesh", &show_mesh);
            if(mesh && mesh->renderer()) {
                mesh->renderer()->set_visible(show_mesh);
            }
            ImGui::SameLine();
            static bool show_sites = false;
            if(sites && sites->renderer()) {
                show_sites = sites->renderer()->is_visible();
            }
            ImGui::Checkbox("sites", &show_sites);
            if(sites && sites->renderer()) {
                sites->renderer()->set_visible(show_sites);
            }
            ImGui::SameLine();
            static bool show_vertices = false;
            if(vertices && vertices->renderer()) {
                show_vertices = vertices->renderer()->is_visible();
            }
            ImGui::Checkbox("vertices", &show_vertices);
            if(vertices && vertices->renderer()) {
                vertices->renderer()->set_visible(show_vertices);
            }
        }
        ImGui::Separator();
        if (ImGui::CollapsingHeader("FilletSeg", ImGuiTreeNodeFlags_DefaultOpen)) {
            float width = 200;
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("eps", &eps, 0.01, 1.0f, "%.2f");
            eps = std::clamp(eps, 0.0, 1.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("radius", &radius, 0.01, 1.0f, "%.2f");
            radius = std::clamp(radius, 0.0, 1.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("s", &s, 0.1, 1.0f, "%.2f");
            s = std::max(s, 1.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("min_score", &min_score, 0.1, 1.0f, "%.2f");
            min_score = std::clamp(min_score, 0.0, 1.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("alpha", &alpha, 0.1, 1.0f, "%.2f");
            alpha = std::clamp(alpha, 0.0, 10.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("angle", &angle, 0.1, 1.0f, "%.2f");
            angle = std::clamp(angle, 0.0, 10.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputInt("nb_neighbors", &nb_neighbors, 1, 1.0f);
            nb_neighbors = std::clamp(nb_neighbors, 0, 100);
            ImGui::SetNextItemWidth(width);
            ImGui::InputInt("num_sor_iter", &num_sor_iter, 1, 1.0f);
            num_sor_iter = std::clamp(num_sor_iter, 0, 5);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("std_ratio", &std_ratio, 0.1, 1.0f, "%.2f");
            std_ratio = std::clamp(std_ratio, 0.0, 1.0);
            ImGui::Separator();
            if (ImGui::Button("scoring", ImVec2(150, 30))) {
                if(mesh) {
                    std::thread t(std::bind(&ViewerImGui::run_scoring, this));
                    t.detach();
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("gcp", ImVec2(150, 30))) {

            }

            if (ImGui::Button("geodesic", ImVec2(310, 30))) {

            }
            if (ImGui::Button("refine_fillet_boundary", ImVec2(310, 30))) {

            }
            if (ImGui::Button("refine_target_normal", ImVec2(310, 30))) {

            }
        }
        ImGui::Separator();
        if (ImGui::CollapsingHeader("DeFillet", ImGuiTreeNodeFlags_DefaultOpen)) {

        }
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    }

    bool ViewerImGui::open() {
        const std::string title("Please choose a file");
        const std::string &default_path =  "../data/";
        const std::vector<std::string> &filters = {
                "Surface Mesh (*.obj *.ply *.off *.stl *.sm *.geojson *.trilist)", "*.obj *.ply *.off *.stl *.sm *.geojson *.trilist"
        };
        const std::vector<std::string> &file_names = dialog::open(title, default_path, filters, true);

        for(auto m : models_) {
            delete_model(m);
        }
        models_.clear();
        int count = 0;
        for (const auto &file_name : file_names) {
            if (add_model(file_name)) {
                mesh = current_model();
                ++count;
                auto now = std::chrono::system_clock::now();
                std::time_t time = std::chrono::system_clock::to_time_t(now);
                std::tm tm;
                localtime_s(&tm, &time); // Windows
                std::ostringstream oss;
                oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
                sim_name = easy3d::file_system::simple_name(file_name);
                base_name = easy3d::file_system::base_name(file_name);

                out_dir = cur_work_dir + "/" + oss.str() +  "_" + base_name + "/";
                scoring_mesh_path = out_dir + "/" + base_name + "_mesh_scoring.ply";
                scoring_sites_path = out_dir + "/" + base_name + "_sites_scoring.ply";
                scoring_vertices_path = out_dir + "/" + base_name + "_vertices_scoring.ply";
                try {
                    // 使用 create_directories 函数创建新的目录（包括父目录）
                    std::filesystem::create_directories(out_dir);
                    std::cout << "Folder created successfully." << std::endl;
                    std::filesystem::copy_file(file_name, out_dir + sim_name);
                } catch (const std::filesystem::filesystem_error& e) {
                    std::cerr << "Failed to create folder: " << e.what() << std::endl;
                }
                break;
            }
        }

        if (count > 0) {
            fit_screen();
            return true;
        }
        return false;
    }

    void ViewerImGui::run_scoring() {
        std::string cli = "scoring.exe -i " + out_dir + sim_name + " -o " + out_dir
                + " -e " + std::to_string(eps) + " -r " + std::to_string(radius)
                + " -s " + std::to_string(s) + " --min_score " + std::to_string(min_score)
                + " --nb_neighbors " + std::to_string(nb_neighbors)
                + " --num_sor_iter " + std::to_string(num_sor_iter)
                + " --std_ratio " + std::to_string(std_ratio);
        std::system(cli.c_str());
        std::cout << "scoring done." <<std::endl;
        state = UPDATE_SCORING;
    }

    void ViewerImGui::update_event() {
        if(state == UPDATE_SCORING) {

        }
    }
}
