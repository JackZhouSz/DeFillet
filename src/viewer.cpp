#include "viewer.h"
#include "my_picker.h"
#include <easy3d/util/dialogs.h>

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
#include <easy3d/renderer/texture_manager.h>
#include <easy3d/gui/picker_surface_mesh.h>
#include <easy3d/renderer/shapes.h>
#include <easy3d/util/timer.h>
#include <easy3d/fileio/resources.h>

#include <3rd_party/imgui/misc/fonts/imgui_fonts_droid_sans.h>
#include <3rd_party/imgui/imgui.h>
#include <3rd_party/imgui/backends/imgui_impl_glfw.h>
#include <3rd_party/imgui/backends/imgui_impl_opengl3.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>


#include <omp.h>

#include <igl/jet.h>
// To have the same shortcut behavior on macOS and other platforms (i.e., Windows, Linux)
#ifdef __APPLE__
#define EASY3D_MOD_CONTROL GLFW_MOD_SUPER
#else
#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#endif



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
            , show_mesh(false)
            , eps(0.03), s(10), radius(0.06)
            , min_score(0.5), alpha(0.5)
            , std_ratio(0.3), num_sor_iter(3), nb_neighbors(30)
            , w_convex(0.08), w_concave(1.0), w1(0.3), w2(0.4)
            , angle(60), beta(1.0), gamma(1.0), num_opt_iter(10), interactive_(false){
        camera()->setUpVector(vec3(0, 1, 0));
        camera()->setViewDirection(vec3(0, 0, -1));
        camera_->showEntireScene();
        fillet_color = easy3d::vec3(1.0, 0, 0);
        non_fillet_color = easy3d::vec3(0.0, 0.0, 1.0);
        selected_color = easy3d::vec3(0.0,1.0,0.0);
        easy3d::vec3 selected_color;
        cur_work_dir = easy3d::file_system::current_working_directory();
        omp_set_num_threads(10);
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


        update_event();


        Viewer::pre_draw();
    }


    void ViewerImGui::post_draw() {
        draw_dashboard();
        draw_corner_axes();
        if (polygon_.size() >= 3) {
            // draw the boundary of the rect/lasso
            shapes::draw_polygon_wire(polygon_, vec4(1.0f, 0.0f, 0.0f, 1.0f), width(), height(), -1.0f);
            // draw its transparent face
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            shapes::draw_polygon_filled(polygon_, vec4(1.0f, 0.0f, 0.0f, 0.2f), width(), height(), -0.9f);
            glDisable(GL_BLEND);
        }
    }


    void ViewerImGui::draw_dashboard(){
        static bool my_tool_active = true;
        static bool tmp = false;
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
        float width = 150;
        if (ImGui::CollapsingHeader("FilletSeg")) {
            if(ImGui::CollapsingHeader("sor_para")) {
                ImGui::SetNextItemWidth(width);
                ImGui::InputInt("nb_neighbors", &nb_neighbors, 1, 1.0f);
                nb_neighbors = std::clamp(nb_neighbors, 0, 100);
                ImGui::SetNextItemWidth(width);
                ImGui::InputInt("num_sor_iter", &num_sor_iter, 1, 1.0f);
                num_sor_iter = std::clamp(num_sor_iter, 0, 5);
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("std_ratio", &std_ratio, 0.1, 1.0f, "%.2f");
                std_ratio = std::clamp(std_ratio, 0.0, 1.0);
            }
            if(ImGui::CollapsingHeader("scoring_para")) {
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
            }
            if(ImGui::CollapsingHeader("gcp_para")) {
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("alpha", &alpha, 0.1, 1.0f, "%.2f");
                alpha = std::clamp(alpha, 0.0, 10.0);
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("w_convex", &w_convex, 0.1, 1.0f, "%.2f");
                w_convex = std::clamp(w_convex, 0.0, 1.0);
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("w_convex", &w_concave, 0.1, 1.0f, "%.2f");
                w_concave = std::clamp(w_concave, 0.0, 1.0);
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("w1", &w1, 0.1, 1.0f, "%.2f");
                w1 = std::clamp(w1, 0.0, 1.0);
                ImGui::SetNextItemWidth(width);
                ImGui::InputDouble("w2", &w2, 0.1, 1.0f, "%.2f");
                w2 = std::clamp(w2, 0.0, 1.0);
            }

            if (ImGui::Button("scoring", ImVec2(150, 30))) {
                if(mesh) {
                    std::thread t(std::bind(&ViewerImGui::run_scoring, this));
                    t.detach();
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("gcp", ImVec2(150, 30))) {
                if(mesh) {
                    std::thread t(std::bind(&ViewerImGui::run_gcp, this));
                    t.detach();
                }
            }


        }
        ImGui::Separator();
        if (ImGui::CollapsingHeader("DeFillet")) {
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("angle", &angle, 0.1, 1.0f, "%.2f");
            angle = std::clamp(angle, 0.0, 180.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("beta", &beta, 1.0, 1.0f, "%.2f");
            beta = std::clamp(beta, 1.0, 10.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputDouble("gamma", &gamma, 1.0, 1.0f, "%.2f");
            gamma = std::clamp(gamma, 1.0, 10.0);
            ImGui::SetNextItemWidth(width);
            ImGui::InputInt("num_opt_iter", &num_opt_iter, 1, 1.0f);
            num_opt_iter = std::clamp(num_opt_iter, 1, 20);
            if (ImGui::Button("geodesic", ImVec2(310, 30))) {
                if(mesh) {
                    std::thread t(std::bind(&ViewerImGui::run_geodesic, this));
                    t.detach();
                }
            }
            if (ImGui::Button("defillet", ImVec2(310, 30))) {
                if(mesh) {
                    std::thread t(std::bind(&ViewerImGui::run_defillet, this));
                    t.detach();
                }
            }
        }
        if (ImGui::CollapsingHeader("Interactive")) {
            if (ImGui::Button("fillet", ImVec2(310, 30))) {
                if(mesh) {
                    update_select(dynamic_cast<easy3d::SurfaceMesh*>(mesh) ,true);
                }
            }
            if (ImGui::Button("non-fillet", ImVec2(310, 30))) {
                if(mesh) {
                    update_select(dynamic_cast<easy3d::SurfaceMesh*>(mesh), false);
                }
            }
            tmp = interactive_;
            ImGui::Checkbox("interactive", &interactive_);
            if(interactive_ != tmp) {
                easy3d::SurfaceMesh* model = dynamic_cast<easy3d::SurfaceMesh*>(mesh);
                auto select = model->face_property<bool>("f:select");
                auto colors = model->face_property<vec3>("f:color");
                auto gcp_label = model->face_property<int>("f:gcp_labels");

                auto drawable = mesh->renderer()->get_triangles_drawable("faces");

                for(auto f : model->faces()) {
                    select[f] = false;
                    colors[f] = gcp_label[f] ? fillet_color : non_fillet_color;
                }
                drawable->set_property_coloring(easy3d::State::FACE, "f:color");
                drawable->update();
            }
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
                scoring_mesh_path = out_dir  + "mesh_scoring.ply";
                scoring_sites_path = out_dir + "sites_scoring.ply";
                scoring_vertices_path = out_dir +  "vertices_scoring.ply";
                gcp_mesh_path = out_dir +  "gcp.ply";
                fillet_geo_path = out_dir + "fillet_geo.ply";
                fiilet_defillet_path = out_dir + "fillet_defillet.ply";
                defillet_path = out_dir + "defillet.ply";
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
        if(std::system(cli.c_str()) == 0) {
            std::cout << "scoring done." << std::endl;
            state = UPDATE_SCORING;
        } else {
            std::cout << "scoring error." << std::endl;
        }
    }

    void ViewerImGui::run_gcp() {

        std::string cli = "gc.exe -i " + scoring_mesh_path + " -o " + out_dir
                          + " -a " + std::to_string(alpha)
                          + " --w_convex " + std::to_string(w_convex)
                          + " --w_concave " + std::to_string(w_concave)
                          + " --w1 " + std::to_string(w1)
                          + " --w2 " + std::to_string(w2);
        if(std::system(cli.c_str()) == 0) {
            std::cout << "gcp done." << std::endl;
            state = UPDATE_GCP;
        } else {
            std::cout << "gcp error" << std::endl;
        }
    }

    void ViewerImGui::run_geodesic() {
        easy3d::io::save_ply(gcp_mesh_path,  dynamic_cast<easy3d::SurfaceMesh*>(mesh), false);
        std::string cli = "geo.exe -i " + gcp_mesh_path + " -o " + out_dir
                          + " --angle " + std::to_string(angle);
        if(std::system(cli.c_str()) == 0) {
            std::cout << "geo done." << std::endl;
            state = UPDATE_GEO;
        } else {
            std::cout << "geo error" << std::endl;
        }
    }

    void ViewerImGui::run_defillet() {
        std::string cli = "defillet.exe -i " + gcp_mesh_path
                          + " -f " + fillet_geo_path
                          + " -o " + out_dir
                          + " --beta " + std::to_string(beta)
                          + " --gamma " + std::to_string(gamma);
                          + " --num_opt_iter " + std::to_string(num_opt_iter);
        if(std::system(cli.c_str()) == 0) {
            std::cout << "defillet done." << std::endl;
            state = UPDATE_DEFILLET;
        } else {
            std::cout << "defillet error" << std::endl;
        }
    }

    void ViewerImGui::update_event() {
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
            model_idx_ = 0;
            state = NOTHING;
        }
        if(state == UPDATE_GCP) {
            for(auto m : models_) {
                if(m) {
                    delete_model(m);
                }
            }
            models_.clear();
            mesh = add_model(gcp_mesh_path);
            sites = add_model(scoring_sites_path);
            vertices = add_model(scoring_vertices_path);
            model_idx_ = 0;
            state = NOTHING;
        }
        if(state == UPDATE_GEO) {
            easy3d::SurfaceMesh* fillet_mesh = easy3d::SurfaceMeshIO::load(fillet_geo_path);
            auto fillet_geo_dis = fillet_mesh->get_vertex_property<float>("v:geo_dis");
            auto fillet_original_index = fillet_mesh->get_vertex_property<int>("v:original_index");
            auto geo_dis = dynamic_cast<easy3d::SurfaceMesh*>(mesh)->vertex_property<float>("v:geo_dis", 0.0);
            for(auto v : fillet_mesh->vertices()) {
                easy3d::SurfaceMesh::Vertex vv(fillet_original_index[v]);
                geo_dis[vv] = fillet_geo_dis[v];
            }
            auto drawable = mesh->renderer()->get_triangles_drawable("faces");
            drawable->set_scalar_coloring(easy3d::State::VERTEX, "v:geo_dis", nullptr, 0.0f, 0.0f);
//            std::cout << EASY3D_RESOURCES_DIR << std::endl;
            const std::string texture_file = SCALAR_COLOR_TEXTURE;
            easy3d::Texture *texture = easy3d::TextureManager::request(texture_file);
            drawable->set_texture(texture);
            drawable->update();

            auto face_tar_normals_x = fillet_mesh->face_property<float>("f:tar_normals_x");
            auto face_tar_normals_y = fillet_mesh->face_property<float>("f:tar_normals_y");
            auto face_tar_normals_z = fillet_mesh->face_property<float>("f:tar_normals_z");
            const easy3d::Box3 &box = fillet_mesh->bounding_box();
            float length = norm(box.max_point() - box.min_point()) * 0.02f;
            std::vector<easy3d::vec3> tmp;
            for(auto f : fillet_mesh->faces()) {
                int num = 0;
                easy3d::vec3 center(0,0,0);
                for(auto v : fillet_mesh->vertices(f)) {
                    center += fillet_mesh->position(v);
                    num++;
                }
                center /= num;
                tmp.emplace_back(center);
                easy3d::vec3 n(face_tar_normals_x[f], face_tar_normals_y[f], face_tar_normals_z[f]);
                tmp.emplace_back(center + n * length);
            }
            easy3d::LinesDrawable* line = mesh->renderer()->add_lines_drawable("target_normals");
            line->update_vertex_buffer(tmp);
            line->set_impostor_type(easy3d::LinesDrawable::ImposterType::CYLINDER);
            line->set_line_width(1.0);
            line->set_uniform_coloring(easy3d::vec4(0.0, 1.0, 0.0, 1.0));
            line->set_visible(true);
            state = NOTHING;
        }
        if(state == UPDATE_DEFILLET) {

            for(auto m : models_) {
                if(m) {
                    delete_model(m);
                }
            }
            models_.clear();
            mesh = add_model(defillet_path);
            sites = add_model(scoring_sites_path);
            vertices = add_model(scoring_vertices_path);
            model_idx_ = 0;
            state = NOTHING;

        }
    }

    bool ViewerImGui::mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) {
        if (interactive_ && mesh) {
            polygon_.push_back(vec2(x, y));
            return false;
        } else
            return Viewer::mouse_drag_event(x, y, dx, dy, button, modifiers);
    }

    bool ViewerImGui::mouse_release_event(int x, int y, int button, int modifiers) {
        if (interactive_ && mesh) {
            if (polygon_.size() >= 3) {

                pick_faces(button,modifiers);
                return true;
            }
            return false;
        } else
            return Viewer::mouse_release_event(x, y, button, modifiers);
    }

    /// Mouse button press event handler
    bool ViewerImGui::mouse_press_event(int x, int y, int button, int modifiers) {
        if (interactive_ && mesh) {
            polygon_.clear();
            polygon_.push_back(vec2(x, y));
            return false;
        } else
            return Viewer::mouse_press_event(x, y, button, modifiers);
    }


    void ViewerImGui::pick_faces(int button, int modifiers) {
        auto mesh_ = dynamic_cast<easy3d::SurfaceMesh*>(mesh);

        if (mesh_) {
            MySurfaceMeshPicker picker(camera());
            picker.pick_faces(mesh_, polygon_, button == GLFW_MOUSE_BUTTON_RIGHT);
            auto select = mesh_->face_property<bool>("f:select");
            auto colors = mesh_->face_property<vec3>("f:color");
            auto gcp_label = mesh_->face_property<int>("f:gcp_labels");


            auto drawable = mesh->renderer()->get_triangles_drawable("faces");

            for(auto f : mesh_->faces()) {
                if(select[f]) {
                    colors[f] = selected_color;
                }
                else {
                    if(gcp_label[f]) {
                        colors[f] = fillet_color;
                    }
                    else {
                        colors[f] = non_fillet_color;
                    }
                }
            }
            drawable->set_property_coloring(easy3d::State::FACE, "f:color");
            drawable->update();

            polygon_.clear();
        }
    }

    void ViewerImGui::update_select(easy3d::SurfaceMesh* model, bool is_fillet) {
        auto select = model->face_property<bool>("f:select");
        auto colors = model->face_property<vec3>("f:color");
        auto gcp_label = model->face_property<int>("f:gcp_labels");

        auto drawable = mesh->renderer()->get_triangles_drawable("faces");

        for(auto f : model->faces()) {
            if(select[f]) {
                gcp_label[f] = is_fillet ? 1 : 0;
                colors[f] = is_fillet ? fillet_color : non_fillet_color;
                select[f] = false;
            }
        }
        drawable->set_property_coloring(easy3d::State::FACE, "f:color");
        drawable->update();
    }
}
