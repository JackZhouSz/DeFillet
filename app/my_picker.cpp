//
// Created by xiaowuga on 3/22/2024.
//

#include "my_picker.h"
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/shader_program.h>
#include <easy3d/renderer/shader_manager.h>
#include <easy3d/renderer/framebuffer_object.h>
#include <easy3d/renderer/opengl_error.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/manipulator.h>
#include <easy3d/util/logging.h>


namespace easy3d {
    MySurfaceMeshPicker::MySurfaceMeshPicker(const Camera *cam) : Picker(cam) {
        use_gpu_if_supported_ = true;
    };

    MySurfaceMeshPicker::~MySurfaceMeshPicker() {
    }

    void MySurfaceMeshPicker::pick_faces(easy3d::SurfaceMesh* model, const Polygon2 &plg, bool deselect) {
        if (!model)
            return;

        ShaderProgram* program = nullptr;
        if (use_gpu_if_supported_) {
            program = ShaderManager::get_program("selection/selection_single_primitive");
            if (!program) {
                std::vector<ShaderProgram::Attribute> attributes;
                attributes.push_back(ShaderProgram::Attribute(ShaderProgram::POSITION, "vtx_position"));
                program = ShaderManager::create_program_from_files("selection/selection_single_primitive", attributes);
            }
            if (!program) {
                use_gpu_if_supported_ = false;
                LOG_N_TIMES(3, ERROR) << "shader program not available, fall back to CPU implementation. " << COUNTER;
            }
        }

        if (use_gpu_if_supported_ && program)
             pick_faces_gpu(model, plg, program, deselect);
    }


    void MySurfaceMeshPicker::pick_faces_gpu(SurfaceMesh *model, const Polygon2 &plg, ShaderProgram* program, bool deselect) {
        auto drawable = model->renderer()->get_triangles_drawable("faces");
        if (!drawable) {
            LOG_N_TIMES(3, WARNING) << "drawable 'faces' does not exist. " << COUNTER;
            return;
        }

        int viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        int width = viewport[2];
        int height = viewport[3];
        setup_framebuffer(width, height);


        const Box2& box = plg.bbox();
        int xmin = box.min_point().x;
        int ymin = box.min_point().y;
        int xmax = box.max_point().x;
        int ymax = box.max_point().y;
        if (xmin > xmax) std::swap(xmin, xmax);
        if (ymin > ymax) std::swap(ymin, ymax);


        std::vector<vec2> region; // the transformed selection region
        for (std::size_t i = 0; i < plg.size(); ++i) {
            const vec2 &p = plg[i];
            float x = p.x;
            float y = p.y;
            region.push_back(vec2(x, y));
        }


        //--------------------------------------------------------------------------
        // render the 'scene' to the new FBO.

        // TODO: the performance can be improved. Since the 'scene' is static, we need to render it to the fbo only
        //       once. Then just query. Re-rendering is needed only when the scene is changed/manipulated, or the size
        //       of the viewer changed.

        // Bind the offscreen fbo for drawing
        fbo_->bind(); easy3d_debug_log_gl_error; easy3d_debug_log_frame_buffer_error;

        float color[4];
        glGetFloatv(GL_COLOR_CLEAR_VALUE, color);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        easy3d_debug_log_gl_error;
        easy3d_debug_log_frame_buffer_error;

        const mat4 MANIP = model->manipulator() ? model->manipulator()->matrix() : mat4::identity();
        program->bind();
        program->set_uniform("MVP", camera()->modelViewProjectionMatrix())
                ->set_uniform("MANIP", MANIP);
        drawable->gl_draw();
        program->release();

        // --- Maybe this is not necessary ---------
        glFlush();
        glFinish();
        // -----------------------------------------


        auto select = model->face_property<bool>("f:select");
//#pragma omp parallel for
        for(int x = xmin; x <= xmax; x++) {
            for (int y = ymin; y <= ymax; y++) {
                int gl_x, gl_y;
                if (geom::point_in_polygon(vec2(x, y), region)) {
                    screen_to_opengl(x, y, gl_x, gl_y, width, height);

                    unsigned char c[4];
                    fbo_->read_color(c, gl_x, gl_y);

                    int id = rgb::rgba(c[0], c[1], c[2], c[3]);
                    if (id >= 0 && id < model->n_faces()) {
                        SurfaceMesh::Face f(id);
                        select[f] = !deselect;
                    }
                }
            }
        }

        fbo_->release(); easy3d_debug_log_gl_error; easy3d_debug_log_frame_buffer_error;
        glClearColor(color[0], color[1], color[2], color[3]);


    }

}