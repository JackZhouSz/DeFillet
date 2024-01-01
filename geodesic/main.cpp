// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
#include "Xin_Wang.h"
#include <iostream>
#include <vector>
using namespace std;
int main(int argc, char* argv[])
{
	CRichModel model(argv[1]);
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
	//for (int i = 0; i < model.GetNumOfVerts(); ++i)
	//	cout << "vert " << i << ": " << alg.GetDistanceField()[i] << endl;
	return 0;
}

