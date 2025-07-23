#include <Inventor/SoDB.h>
#include <Inventor/SoInput.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoMFInt32.h>
#include <Inventor/fields/SoMFColor.h>
#include <iostream>

void extractMeshData(SoNode* node) {
    // 确保节点是 SoIndexedFaceSet 类型
    if (node->isOfType(SoIndexedFaceSet::getClassTypeId())) {
        SoIndexedFaceSet* ifs = (SoIndexedFaceSet*)node;
        // return;
        // 获取顶点数据
        SoVertexProperty* vertexProperty = ifs->vertexProperty.getValue();
        auto asd = vertices.getValue(0);
        // std::cout << "Vertices (" << vertices.getNum() << "):" << std::endl;
        // for (int i = 0; i < vertices.getNum(); ++i) {
        //     const SbVec3f& vertex = vertices[i];
        //     std::cout << "  " << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << std::endl;
        // }

        // 获取颜色数据（如果有的话）
        // SoMFColor& colors = ifs->vertexProperty->orderedRGBA;
        // std::cout << "Colors (" << colors.getNum() << "):" << std::endl;
        // for (int i = 0; i < colors.getNum(); ++i) {
        //     const SbColor& color = colors[i];
        //     std::cout << "  " << color[0] << ", " << color[1] << ", " << color[2] << std::endl;
        // }

        // 获取面的索引数据
        SoMFInt32& coordIndex = ifs->coordIndex;
        std::cout << "Face Indices (" << coordIndex.getNum() << "):" << std::endl;
        for (int i = 0; i < coordIndex.getNum(); ++i) {
            int index = coordIndex[i];
            if (index == -1) {
                std::cout << "  End of face" << std::endl;  // -1 表示面结束
            } else {
                std::cout << "  " << index << std::endl;
            }
        }

        // 获取材质索引数据
        SoMFInt32& materialIndex = ifs->materialIndex;
        std::cout << "Material Indices (" << materialIndex.getNum() << "):" << std::endl;
        for (int i = 0; i < materialIndex.getNum(); ++i) {
            std::cout << "  " << materialIndex[i] << std::endl;
        }
    }
}

void traverseScene(SoNode* node) {
    // 如果是一个容器节点（如 SoGroup、SoSeparator），递归遍历所有子节点
    if (node->isOfType(SoGroup::getClassTypeId())) {
        SoGroup* group = (SoGroup*)node;
        for (int i = 0; i < group->getNumChildren(); ++i) {
            traverseScene(group->getChild(i));
        }
    } else {
        // 如果是 SoIndexedFaceSet 类型，则提取网格数据
        extractMeshData(node);
    }
}

int main() {
    // 初始化 Coin3D 库
    SoDB::init();

    // 创建输入对象
    SoInput input;

    // 打开指定的 .iv 文件
    if (!input.openFile("E:\\DeFillet\\Fillet Identification\\rgb\\HFP132.iv")) {
        std::cerr << "Unable to open .iv file!" << std::endl;
        return 1;
    }

    // 读取文件内容
    SoSeparator* root = SoDB::readAll(&input);
    if (root == nullptr) {
        std::cerr << "Failed to parse .iv file!" << std::endl;
        return 1;
    }

    // 遍历根节点并提取网格数据
    traverseScene(root);

    // 清理内存
    SoDB::finish();

    return 0;
}

// "E:\\DeFillet\\Fillet Identification\\rgb\\HFP132.iv"