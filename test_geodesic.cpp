//
// Created by xiaowuga on 12/27/2023.
//
//#include "Xin_Wang.h"
#include <Xin_Wang.h>
#include <iostream>
#include <vector>
using namespace std;
int main()
{
    CRichModel model("D:\\code\\defillet\\data\\bottle_fillet_remeshed.obj");
    model.LoadModel();

    int id1 = 0;
    int id2 = 100;
    int id3 = 5;
    set<int> sources;
    sources.insert(id1);
    sources.insert(id2);
    sources.insert(id3);

    CXin_Wang alg(model, sources);
    alg.Execute();

    //observe distances...
//    for (int i = 0; i < model.GetNumOfVerts(); ++i) {
//        cout << "vert " << i << ": " << alg.GetDistanceField()[i] << endl;
//    }
    for (int i = 0; i < 10; ++i) {
        cout << "vert " << i << ": " << alg.GetDistanceField()[i] << endl;
        cout << alg.GetAncestor(i) << std::endl;
    }
    return 0;
}