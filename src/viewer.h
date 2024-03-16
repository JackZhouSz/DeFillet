#ifndef EASY3D_TUTORIAL_IMGUI_VIEWER_H
#define EASY3D_TUTORIAL_IMGUI_VIEWER_H


#include <easy3d/viewer/viewer.h>

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <fillet_seg.h>
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
            NOTHING
        };

    protected:

        // imgui plugins
        void init() override;

        // draw the widgets
        void pre_draw() override;

        void draw() const override;

        //  the widgets
        void post_draw() override;

        void cleanup() override;

        void post_resize(int w, int h) override;

        void update_event();

        bool callback_event_cursor_pos(double x, double y) override;
        bool callback_event_mouse_button(int button, int action, int modifiers) override;
        bool callback_event_keyboard(int key, int action, int modifiers) override;
        bool callback_event_character(unsigned int codepoint) override;
        bool callback_event_scroll(double dx, double dy) override;

        void draw_dashboard();


        void run_scoring();
        void run_gcp();
        void run_geodesic();
        void run_refine_fillet_boundary();
        void run_refine_target_normal();

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
        FilletSeg* fillet_seg;
        easy3d::SurfaceMesh* mesh;
        easy3d::PointCloud* sites;
        easy3d::PointCloud* vertices;
        bool show_mesh;
        double eps;
        double radius;
        double s;
        double min_score;
        double alpha;
        double angle;
        STATE state;
    };

}

#endif	// EASY3D_TUTORIAL_IMGUI_VIEWER_H
