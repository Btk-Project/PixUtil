#include <SDL2/SDL.h>

#define PIXUTIL_SDL_EXTERNAL

// #define PIXUTIL_MALLOC SDL_malloc
// #define PIXUTIL_REALLOC SDL_realloc
// #define PIXUTIL_FREE SDL_free


#include "../pixutil_raster.hpp"
#include "../pixutil_filter.hpp"
#include "../pixutil.hpp"
using namespace PixRaster;
using namespace PixFilter;
using namespace PixUtil;

int main(){
    // Create a window and renderer
    SDL_Window *window = SDL_CreateWindow("Test",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        640,480,
        SDL_WINDOW_SHOWN
    );
    // SDL_Renderer *renderer = SDL_CreateRenderer(window,-1,SDL_RENDERER_ACCELERATED);
    // SDL_Texture *texture = SDL_CreateTexture(renderer,SDL_PIXELFORMAT_RGBA32,SDL_TEXTUREACCESS_STREAMING,640,480);    

    // SDL_Event event;
    // int x = 0,y = 0;
    // while(SDL_WaitEvent(&event)){
    //     if(event.type == SDL_QUIT){
    //         break;
    //     }
    //     //Draw a line by mouse
    //     if(event.type == SDL_MOUSEMOTION){
    //         x = event.motion.x;
    //         y = event.motion.y;
    //     }

    //     SDL_SetRenderDrawColor(renderer,0,0,0,255);
    //     SDL_RenderClear(renderer);
    //     //Shade at texture
    //     void *pixels;
    //     int pitch;
    //     SDL_LockTexture(texture,nullptr,&pixels,&pitch);
    //     memset(pixels,0,640 * 480 * 4);

    //     RGBAView view(pixels,640,480,pitch);

    //     auto shade = AddScissor(100,100,540,380,[&](int x,int y){
    //         view.at(x,y) = Color{255,0,0,255};
    //     });

    //     DrawLine(0,0,x,y,shade);
    //     DrawCircle(x,y,50,shade);
    //     SDL_UnlockTexture(texture);


    //     SDL_RenderCopy(renderer,texture,nullptr,nullptr);
    //     SDL_RenderPresent(renderer);
    // }

    SDL_Surface *surface = SDL_GetWindowSurface(window);
    SDLSurfaceView view(surface);

    SDL_Event event;
    int x = 0,y = 0;
    while(SDL_WaitEvent(&event)){
        if(event.type == SDL_QUIT){
            break;
        }
        if(event.type == SDL_MOUSEMOTION){
            x = event.motion.x;
            y = event.motion.y;
        }

        FillRect(view,0,0,640,480,Color{0,0,0,255});

        auto shade = AddScissor(0,0,640,480,[&](int x,int y){
            view.at(x,y) = Color{
                Uint8(rand() % 255),
                Uint8(rand() % 255),
                Uint8(rand() % 255),
                255
            };
        });

        DrawLine(0,0,x,y,shade);

        DrawCircle(x,y,50,shade);
        //Top left
        
        FillTriangle(x+100,y + 100,x - 100,y - 200,x - 60,y + 200,shade);

        SDL_UpdateWindowSurface(window);
    }

    //Cleanup
    // SDL_DestroyTexture(texture);
    // SDL_DestroyRenderer(renderer);
    SDL_FreeSurface(surface);
    SDL_DestroyWindow(window);
    SDL_Quit();
}