#define PIXUTIL_SDL_EXTERNAL
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <iostream>
#define PIXUTIL_MATRIX_IO
#include "../pixutil_filter.hpp"
#include "../pixutil.hpp"

int main(int argc,char **argv){
    SDL_Surface *surf = IMG_Load(
        argv[1]
    );
    if(surf == nullptr){
        printf("Load failed\n");
        return -1;
    }
    SDL_Surface *out = SDL_DuplicateSurface(surf);
    //Gray buffer
    PixUtil::SDLSurfaceView view(surf);
    PixUtil::SDLSurfaceView dst(out);

    PixFilter::DFT(dst,view);

    PixUtil::SDLShowView(dst);
    PixUtil::SDLShowView(view);

    PixUtil::SDLWaitForQuit();

    SDL_FreeSurface(surf);
    SDL_FreeSurface(out);
    SDL_Quit();
}