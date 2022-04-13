#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <iostream>

#define PIXUTIL_MATRIX_IO
#define PIXUTIL_SDL_EXTERNAL
#include "../pixutil_filter.hpp"
#include "../pixutil.hpp"

using namespace PixUtil;
using namespace PixFilter;

int main(int argc,char **argv){
    SDL_Surface *src_surface = IMG_Load(argv[1]);
    SDL_Surface *dst_surface = SDL_CreateRGBSurfaceWithFormat(
        0,
        100,
        100,
        32,
        SDL_PIXELFORMAT_ARGB8888
    );
    SDL_Surface *rot_surface = SDL_DuplicateSurface(src_surface);

    SDLSurfaceView src_view(src_surface);
    SDLSurfaceView dst_view(dst_surface);
    SDLSurfaceView rot_view(rot_surface);

    BilinearScale(dst_view,src_view);


    rot_view.fill(0);
    Rotate(rot_view,src_view,290);

    //Show
    SDLShowView(src_view);
    SDLShowView(dst_view);
    SDLShowView(rot_view);

    //Wait
    SDLWaitForQuit();
    SDL_SaveBMP(rot_surface,"big_out.bmp");

    //Free
    SDL_FreeSurface(src_surface);
    SDL_FreeSurface(dst_surface);
    SDL_FreeSurface(rot_surface);
    SDL_Quit();
}