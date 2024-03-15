#include <easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/random.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
//#include <easy3d/util/initializer.h>
#include <easy3d/util/timer.h>


using namespace easy3d;

// This example shows how to use another thread for
//      - repeatedly modifying a model, and
//      - notifying the viewer thread


// a function that modifies the model.
// in this simple example, we add more points (with per point colors) to a point cloud.
void edit_model(PointCloud *cloud, Viewer *viewer) {
    if (cloud->n_vertices() >= 1000000) // stop growing when the model is too big
        return;

    auto colors = cloud->vertex_property<vec3>("v:color");
    for (int i = 0; i < 100; ++i) {
        auto v = cloud->add_vertex(vec3(random_float(), random_float(), random_float()));
        colors[v] = vec3(random_float(), random_float(), random_float()); // we use a random color
    }

    // notify the renderer to update the OpenGL buffers
    cloud->renderer()->update();
    // notify the viewer to update the display
    viewer->update();

    std::cout << "#points: " << cloud->n_vertices() << std::endl;
}


int main(int argc, char **argv) {
    // initialize Easy3D.
//    initialize();

    // create the viewer
    Viewer viewer("ASDF");

    // create a point cloud (with a per point color property) from a set of random points
    auto cloud = new PointCloud;
    // add a per point color property
    auto colors = cloud->add_vertex_property<vec3>("v:color");
    for (int i = 0; i < 100; ++i) {
        auto v = cloud->add_vertex(easy3d::vec3(random_float(), random_float(), random_float()));
        // assign a color to each point (here simply a red color)
        colors[v] = vec3(1.0f, 0.0f, 0.0f);
    }
    // add the point cloud to the viewer for visualization
    viewer.add_model(cloud);

    // set up visualization parameters
    auto drawable = cloud->renderer()->get_points_drawable("vertices");
    // set point size
    drawable->set_point_size(10.0f);
    // set coloring method: we want to visualize the point cloud using the per point color property
    drawable->set_coloring(Drawable::COLOR_PROPERTY, Drawable::VERTEX, "v:color");

    // run the process in another thread.
#if 1   //  we use a timer to repeatedly edit the point cloud every 300 milliseconds
    Timer<PointCloud *, Viewer *> timer;
    timer.set_interval(300, edit_model, cloud, &viewer);

    // stop editing the model after 20 seconds
    Timer<>::single_shot(20000, &timer, &Timer<PointCloud*, Viewer*>::stop); // or Timer<>::single_shot(4900, [&]() -> void { timer.stop(); });
#else   // use a simple for loop to repeat the editing a fix number of iterations
    Timer<>::single_shot(0, [&]() {
        // in this example, we modify the point cloud 50 times
        for (int i = 0; i < 50; ++i) {
            // allow sufficient time for each edit to complete (each edit runs in a separate thread)
            std::this_thread::sleep_for(std::chrono::milliseconds(300));
            std::thread t(edit_model, cloud, &viewer);
            t.detach();
        }
    });
#endif

    // run the viewer
    return viewer.run();
}
