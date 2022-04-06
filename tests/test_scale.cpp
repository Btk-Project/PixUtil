#include <SDL2/SDL.h>

#define PIXUTIL_SDL_EXTERNAL
#include "../pixutil_filter.hpp"
#include "../pixutil.hpp"

using namespace PixUtil;
using namespace PixFilter;

int main(){
    SDL_Surface *src_surface = SDL_LoadBMP("test.bmp");
    SDL_Surface *dst_surface = SDL_CreateRGBSurfaceWithFormat(
        0,
        100,
        100,
        32,
        SDL_PIXELFORMAT_ARGB8888
    );

    SDLSurfaceView src_view(src_surface);
    SDLSurfaceView dst_view(dst_surface);

    BilinearScale(dst_view,src_view);
    GaussianBlur(dst_view,2,4);

    //Show
    // SDLShowView(src_view);
    SDLShowView(dst_view);

    //Wait
    SDLWaitForQuit();
    // SDL_SaveBMP(dst_surface,"big_out.bmp");

    //Free
    SDL_FreeSurface(src_surface);
    SDL_FreeSurface(dst_surface);
    SDL_Quit();
}