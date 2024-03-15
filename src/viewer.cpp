#include "viewer.h"
#include "easy3d/util/dialogs.h"

#include <iostream>
#include <cstdio>
#include <thread>

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
            , eps(0.03), s(10), radius(0.1)
            , min_score(0.5), alpha(0.5) {
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
        ImGui::SetNextWindowSize(ImVec2(325, h), ImGuiCond_Always);
        ImGui::SetNextWindowPos(ImVec2(w - 325, 0), ImGuiCond_Always);
        ImGui::Begin("Dashboard", &my_tool_active, ImGuiWindowFlags_MenuBar);
        ImGui::SetNextWindowContentSize(ImVec2(0, 0));
        if (ImGui::CollapsingHeader("Open/Save", ImGuiTreeNodeFlags_DefaultOpen)) {

            if (ImGui::Button("Open", ImVec2(150, 30))) {
                const std::string title("Please choose a file");
                const std::string &default_path = "../data/";
                const std::vector<std::string> &filters = {
                        "Surface Mesh (*.obj *.ply)", "*.obj *.ply"
                };
                const std::vector<std::string> &file_names = easy3d::dialog::open(title, default_path, filters, true);
                if(file_names.size() > 0) {
                    if (mesh) {
                        mesh->renderer()->set_visible(false);
                        fillet_seg->reset();
                        delete mesh->renderer();
                        delete mesh->manipulator();
                        delete mesh;
                    }
                    if(fillet_seg) {
                        delete fillet_seg;
                    }
                    fillet_seg = new FilletSeg();
                    mesh = easy3d::SurfaceMeshIO::load(file_names[0]);
                    auto renderer = new easy3d::Renderer(mesh, true);
                    mesh->set_renderer(renderer);
                    auto manipulator = new Manipulator(mesh);
                    mesh->set_manipulator(manipulator);
                    fillet_seg->set_mesh(mesh);
                    show_mesh = true;
                    if (sites) {
                        sites->renderer()->set_visible(false);
                        delete sites->renderer();
                        delete sites->manipulator();
                        delete sites;
                    }
                    if (vertices) {
                        vertices->renderer()->set_visible(false);
                        delete vertices->renderer();
                        delete vertices->manipulator();
                        delete vertices;
                    }
                    Viewer::fit_screen(mesh);
                }
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

            ImGui::Separator();
            if (ImGui::Button("scoring", ImVec2(150, 30))) {
                if(mesh && fillet_seg) {
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

    void ViewerImGui::run_scoring() {

        fillet_seg->reset();
        fillet_seg->set_mesh(mesh);
        fillet_seg->set_eps(eps);fillet_seg->set_radius(radius);
        fillet_seg->set_min_score(min_score);
        fillet_seg->set_s(s);
        fillet_seg->run_scoring();
        std::cout << "ASD" <<std::endl;
        mesh->renderer()->set_visible(false);
        auto scores = mesh->face_property<float>("f:scores");
        int nb_faces = mesh->n_faces();
        Eigen::VectorXd Z(nb_faces);
        for (int i = 0; i < nb_faces; i++) {
            Z[i] = scores.vector()[i];
        }
        Eigen::MatrixXd Ct;
        igl::jet(Z, true, Ct);
        auto coloring_mesh = mesh->face_property<vec3>("f:color", vec3(0, 0, 0));
        for(auto f : mesh->faces()) {
            coloring_mesh[f] = easy3d::vec3(Ct(f.idx(), 0),
                                       Ct(f.idx(), 1), Ct(f.idx(), 2));
        }
        auto drawable_mesh = mesh->renderer()->get_triangles_drawable("faces");
        drawable_mesh->set_property_coloring(State::FACE, "f:color");

        if (sites) {
            delete sites->renderer();
            delete sites->manipulator();
            delete sites;
        }
        sites = new easy3d::PointCloud;
        auto renderer_sites = new easy3d::Renderer(sites, true);
        sites->set_renderer(renderer_sites);
        auto manipulator_sites = new Manipulator(sites);
        sites->set_manipulator(manipulator_sites);
        sites->renderer()->set_visible(false);
        auto& tmp_sites = fillet_seg->get_sites();
        for(auto v : tmp_sites) {
            sites->add_vertex(v);
        }
        auto coloring_sites= sites->vertex_property<vec3>("v:color", vec3(0, 0, 0));
        for(auto v : sites->vertices()) {
            coloring_sites[v] = easy3d::vec3(Ct(v.idx(), 0),
                                            Ct(v.idx(), 1), Ct(v.idx(), 2));
        }
        auto drawable_sites = sites->renderer()->get_points_drawable("vertices");
        drawable_sites->set_property_coloring(State::VERTEX, "v:color");


        if (vertices) {
            delete vertices->renderer();
            delete vertices->manipulator();
            delete vertices;
        }
        vertices = new easy3d::PointCloud;
        auto renderer_vertices = new easy3d::Renderer(vertices, true);
        vertices->set_renderer(renderer_vertices);
        auto manipulator_vertices = new Manipulator(vertices);
        vertices->set_manipulator(manipulator_vertices);
        vertices->renderer()->set_visible(false);
        auto& tmp_vertices = fillet_seg->get_vertices();
        auto& tmp_vertices_scores = fillet_seg->get_vertices_scores();
        auto& tmp_index = fillet_seg->get_vertices_index_inside_box();

        int num = tmp_index.size();
        Z.resize(num);
        for(int i = 0; i < num; i++) {
            vertices->add_vertex(tmp_vertices[tmp_index[i]]);
            Z[i] = tmp_vertices_scores[tmp_index[i]];
        }
        igl::jet(Z, true, Ct);
        auto drawable_vertices = vertices->renderer()->get_points_drawable("vertices");
        drawable_vertices->set_property_coloring(State::VERTEX, "v:color");
        auto coloring_vertices = vertices->vertex_property<vec3>("v:color", vec3(0, 0, 0));
        for(auto v : vertices->vertices()) {
            coloring_vertices[v] = easy3d::vec3(Ct(v.idx(), 0),
                                                Ct(v.idx(), 1), Ct(v.idx(), 2));
        }
        drawable_vertices->set_property_coloring(State::VERTEX, "v:color");


        drawable_sites->update();
        drawable_mesh->update();
        drawable_vertices->update();
        mesh->renderer()->set_visible(true);
        update();
    }
}
