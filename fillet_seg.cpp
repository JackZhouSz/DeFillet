//
// Created by xiaowuga on 3/11/2024.
//

#include "fillet_seg.h"
#include "MeshVoronoi3d.h"

void FilletSeg::seg() {
    site_scoring();
    run_gc();
}

void