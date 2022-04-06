#include "../pixutil_image.hpp"

int main(int argc,char **argv){
    auto bitmap = PixImage::LoadFromFilename(argv[1]);
    if(bitmap.empty()){
        printf("Load failed\n");
        return -1;
    }
    PIXUTIL_FREE(bitmap.pixels);
}