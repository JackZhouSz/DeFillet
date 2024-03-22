//
// Created by xiaowuga on 3/22/2024.
//

#ifndef DEFILLET_MY_PICKER_H
#define DEFILLET_MY_PICKER_H


#include <easy3d/gui/picker.h>
#include <easy3d/core/surface_mesh.h>

namespace easy3d {
    class ShaderProgram;

    class MySurfaceMeshPicker : Picker {
    public:
        MySurfaceMeshPicker(const Camera* cam);
        ~MySurfaceMeshPicker();

        void pick_faces(easy3d::SurfaceMesh* model, const Polygon2 &plg, bool deselect);

        void pick_faces_gpu(easy3d::SurfaceMesh* model, const Polygon2 &plg, ShaderProgram* program, bool deselect);

    private:
        bool use_gpu_if_supported_;

    };
}


#endif //DEFILLET_MY_PICKER_H
