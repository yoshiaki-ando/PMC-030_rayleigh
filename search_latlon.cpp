#include <iostream>
#include <fstream>
#include <sstream>

double search_latlon(
  double lat,
  int side
  ){

  std::ifstream ifs("Latitude_Longitude.dat");//data.datを開く

  //ファイルが開けない場合
  if(ifs.fail()){
    std::cout << "入力ファイルをオープンできませんでした" << std::endl;
    //return 1;
  }

  std::string line,tmp;
  int i=0;
  double lat_if=0.0,lon=0;
  while (getline(ifs,line)) {
    std::istringstream stream(line);
    while(getline(stream,tmp,' ')){
      if(i>1 && i%2==0){
        lat_if = std::stod(tmp);
        if(std::abs(lat_if - lat) < 0.1){
          getline(stream,tmp,' ');
          lon = std::stod(tmp);
          if (side > 0) {
            break;//東(140°E 以下)ならbreakあり,西(140°E以上)ならbreakなし
          }
        }
      }
    i++;
    }
  }
  return lon;
}
