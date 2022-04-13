#include <ft2build.h>
#include <string>

#include FT_FREETYPE_H

#include "../pixutil.hpp"

static FT_Library library;
static FT_Face face;


auto measure(const std::u32string &str) -> std::pair<int,int>{
    int w = 0;
    int h = 0;
    for(auto c : str){
        FT_Load_Char(face,c,FT_LOAD_RENDER);
        w += face->glyph->bitmap.width;
        h = std::max<int>(h,face->glyph->bitmap.rows);
    }
    return std::make_pair(w,h);
}

int main(){
    FT_Init_FreeType(&library);
    FT_New_Face(library,"/usr/share/fonts/truetype/noto/NotoSansLaoUI-ExtraCondensedExtraLight.ttf",0,&face);

    auto s = measure(U"HelloWorld");

    FT_Done_FreeType(library);
}