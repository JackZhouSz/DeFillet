//
// Created by xiaowuga on 2024/11/22.
//

#ifndef VORO_VIEWER_H
#define VORO_VIEWER_H

#include <easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>

namespace easy3d{

    class VoroViewer : public Viewer {
    public:
        VoroViewer(const std::string& title = "VoroViewer") : Viewer(title){

        }
        void init(std::vector<easy3d::vec3>& vor_points
          , easy3d::SurfaceMesh* mesh
          , std::vector<std::vector<int>>& vvns
          , std::vector<std::vector<int>>& vvcs
          , std::vector<std::vector<int>>& scvv);
        void mark_selection1(PointCloud::Vertex v);
        void mark_selection2(SurfaceMesh::Face f);
    private:
        /// Mouse button press event handler
        bool mouse_press_event(int x, int y, int button, int modifiers) override;
    private:
        easy3d::PointCloud* cloud_;
        easy3d::SurfaceMesh* mesh_;
        std::vector<std::vector<int>> vvns_;
        std::vector<std::vector<int>> vvcs_;
        std::vector<std::vector<int>> scvv_;
    };
}



#endif //VORO_VIEWER_H
