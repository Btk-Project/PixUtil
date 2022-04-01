#include "../pixutil_raster.hpp"
#include "../pixutil.hpp"

using namespace PixUtil;
using namespace PixRaster;

int main(){
    //Random generate a bitmap
    int w = 10;
    int h = 10;
    int stride = w * 4;
    Uint8 *pixels = new Uint8[stride * h];
    for(int i = 0;i < h;i++){
        for(int j = 0;j < w;j++){
            pixels[i * stride + j * 4 + 0] = rand() % 256;
            pixels[i * stride + j * 4 + 1] = rand() % 256;
            pixels[i * stride + j * 4 + 2] = rand() % 256;
            pixels[i * stride + j * 4 + 3] = rand() % 256;
        }
    }
    //Create a view
    RGBAView view(pixels,w,h);
    //Test locating
    for(auto row:view){
        for(auto pix:row){
            assert(
                pix == reinterpret_cast<uint32_t*>(pixels)[pix.y * w + pix.x]
            );
        }
        printf("\n");
    }
    FillRect(view,0,1,100,10,0);
    view[0][0] = Color{255,0,0,255};
    for(auto row:view){
        for(auto pix:row){
            printf("%u ",Uint32(pix));
        }
        printf("\n");
    }

    //Test drawing
    auto shade = [&](int x,int y){
        view.at(x,y) = Color{255,0,0,255};
    };
    DrawLine(0,0,10,10,AddScissor(0,0,10,10,shade));

    //Cleanup
    delete[] pixels;
}