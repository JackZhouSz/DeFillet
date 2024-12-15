//
// Created by xiaowuga on 2024/11/22.
//

#ifndef VORO_VIEWER_H
#define VORO_VIEWER_H

#include <easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>

namespace easy3d{
    class VoroViewer : public Viewer {
    public:
        VoroViewer(const std::string& title = "VoroViewer") : Viewer(title){
            pd = new PointsDrawable("pd");
            add_drawable(pd);
            pd->set_visible(false);
        }
        void init(std::vector<easy3d::vec3>& vor_points
          , easy3d::SurfaceMesh* mesh
          , std::vector<std::vector<int>>& vvns
          , std::vector<std::vector<int>>& vvcs);

        void mark_selection1(PointCloud::Vertex v);
    private:
        /// Mouse button press event handler
        bool mouse_press_event(int x, int y, int button, int modifiers) override;
        bool key_press_event(int key, int modifiers) override;
    private:
        easy3d::PointCloud* cloud_;
        easy3d::SurfaceMesh* mesh_;
        std::vector<std::vector<int>> vvns_;
        std::vector<std::vector<int>> vvcs_;
        std::vector<std::set<int>> scvv_;
        PointsDrawable* pd;
    };
}



#endif //VORO_VIEWER_H
