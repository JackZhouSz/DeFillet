#ifndef EASY3D_TUTORIAL_IMGUI_VIEWER_H1
#define EASY3D_TUTORIAL_IMGUI_VIEWER_H1


#include <easy3d/viewer/viewer.h>

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
//#dkm "./fillet_seg.h"
// A very good tutorial for imgui:
// https://eliasdaler.github.io/using-imgui-with-sfml-pt1/
// https://eliasdaler.github.io/using-imgui-with-sfml-pt2/


struct ImGuiContext;


namespace easy3d {

    class ViewerImGui : public Viewer
    {
    public:
        ViewerImGui(
                const std::string& title = "Easy3D ImGui Viewer",
                int samples = 4,
                int gl_major = 3,
                int gl_minor = 2,
                bool full_screen = false,
                bool resizable = true,
                int depth_bits = 24,
                int stencil_bits = 8
        );

        enum STATE{
            UPDATE_SCORING,
            UPDATE_GCP,
            UPDATE_GEO,
            UPDATE_DEFILLET,
            NOTHING
        };

    protected:

        // imgui plugins
        void init() override;

        // draw the widgets
        void pre_draw() override;

        //  the widgets
        void post_draw() override;

        void cleanup() override;

        void post_resize(int w, int h) override;

        bool open() override;

        void update_event();

        void pick_faces(int button, int modifiers);

        void update_select(easy3d::SurfaceMesh* model, bool is_fillet);

        bool callback_event_cursor_pos(double x, double y) override;
        bool callback_event_mouse_button(int button, int action, int modifiers) override;
        bool callback_event_keyboard(int key, int action, int modifiers) override;
        bool callback_event_character(unsigned int codepoint) override;
        bool callback_event_scroll(double dx, double dy) override;
        /// Mouse button press event handler
        bool mouse_press_event(int x, int y, int button, int modifiers) override;
        /// Mouse button release event handler
        bool mouse_release_event(int x, int y, int button, int modifiers) override;
        /// Mouse drag (i.e., a mouse button was pressed) event handler
        bool mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) override;


        void draw_dashboard();


        void run_scoring();
        void run_gcp();
        void run_geodesic();
        void run_defillet();

    protected:
        // Ratio between the framebuffer size and the window size.
        // May be different from the DPI scaling!
        double pixel_ratio();

        double widget_scaling() { return dpi_scaling() / pixel_ratio(); }

        // We don't need a per-window font. So this function is static
        void  reload_font(int font_size = 16);

    protected:
        // Single global context by default, but can be overridden by the user
        static ImGuiContext *	context_;

        // Global variables for all the windows
        float	alpha_;
        bool	movable_;

        float   menu_height_;


    private:
        std::string cur_work_dir;
        std::string out_dir;
        std::string sim_name;
        std::string base_name;
        std::string scoring_mesh_path;
        std::string scoring_sites_path;
        std::string scoring_vertices_path;
        std::string gcp_mesh_path;
        std::string fillet_geo_path;
        std::string fiilet_defillet_path;
        std::string defillet_path;
        easy3d::Model* mesh;
        easy3d::Model* sites;
        easy3d::Model* vertices;
        bool show_mesh;
        double eps;
        double radius;
        double s;
        double min_score;
        double alpha;
        int nb_neighbors;
        int num_sor_iter;
        double std_ratio;
        double w_convex;
        double w_concave;
        double w1, w2;

        double angle;
        double beta;
        double gamma;
        int num_opt_iter;
        STATE state;

        easy3d::Polygon2 polygon_;
        bool interactive_;
        easy3d::vec3 fillet_color;
        easy3d::vec3 non_fillet_color;
        easy3d::vec3 selected_color;
    };

}

#endif	// EASY3D_TUTORIAL_IMGUI_VIEWER_H
