#ifndef __Map_h__
#define __Map_h__

#include <vector>
#include "random.h"

struct Pos{
  double x;
  double y;
};

class Map{

public:
  Map();
	Map(int, int, Random*, double);
	~Map();

  double GetX(int);
  double GetY(int);
  double sq_distance(int,int);

protected:
	int _Nc;							//number of cities
	std::vector<Pos> _cities;
};

#endif
