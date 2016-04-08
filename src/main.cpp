#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>
#include <cstdio>
#include <list>

#define MyINFINITY 4294967295

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;
using std::vector;
using std::iterator;
using std::list;
using std::tuple;

#include "io.h"
#include "matrix.h"
#include "MyObject.h"


class SimpleDT{
private:
	double *ans = NULL;
	uint size = 0;
public:
	SimpleDT(const vector<double>& input)
	{
		ans = new double[input.size()];
		size = input.size();
		implamatation_of_distance_transform(input);
	}
	void implamatation_of_distance_transform(const vector<double> &f){
		uint k = 0;
		double *v = new double[f.size()], *z = new double[f.size() + 10];
		v[0] = 0;
		z[0] = -MyINFINITY;
		z[1] = MyINFINITY;
		double s;
		for (uint q = 1; q < f.size(); q++){
			while ((s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k]))/(2 * q - 2 * v[k])) <= z[k]){
				k--;
			}
			k++;
			v[k] = q;
			z[k] = s;
			z[k + 1] = MyINFINITY;
		}
		k = 0;
		for (uint q = 0; q < f.size(); q++){
			while (z[k + 1] < q)
				k++;
			ans[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
		}
		delete []z;
		delete []v;					
	}
	
	~SimpleDT (){
		delete []ans;
	}
	
	SimpleDT(const SimpleDT& o){
		ans = new double[o.size];
		size = o.size;
		memcpy(ans, o.ans, sizeof(double) * o.size);
	}
	
	SimpleDT &operator= (const SimpleDT &o){
		ans = new double[o.size];
		size = o.size;
		memcpy(ans, o.ans, sizeof(double) * o.size);
		return *this;
	}
	
	double distance_transform(uint q){
		return ans[q];
	}
};
		
		
vector<double> Matrix2vector(const Matrix<uint>& in, uint q){
	vector<double> ans;
	for (uint i = 0; i < in.n_rows; i++){
		ans.push_back( double(in(i, q)));
	}
	return ans;
}


class DistanceTransform
{
	Matrix<uint> grind;
	Matrix<uint> answer;
public:
	DistanceTransform (const Matrix<uint> &_in): grind(_in.deep_copy()), answer(_in){
		for (uint i = 0; i < grind.n_rows; i++){
			for (uint j = 0; j < grind.n_cols; j++){
				if (grind(i,j) != 0){
					grind(i,j) = MyINFINITY;
				}
			}
		}
	}
	void distance_transform(){
		vector<SimpleDT> hlp;
		for (uint i = 0; i < grind.n_cols; i++)
		{
			hlp.push_back(SimpleDT(Matrix2vector(grind, i)));
		}
		vector<double> tmp(grind.n_cols, 0);
	
		for (uint i = 0; i < grind.n_rows; i++){
			for(uint j = 0; j < grind.n_cols; j++){
				tmp[j] = hlp[j].distance_transform(i);
			}
			SimpleDT temp2DT(tmp);
			for (uint j = 0; j < grind.n_cols; j++){
				answer(i,j) = temp2DT.distance_transform(j);
			}
		}
	}
	Matrix<uint> &get_answer(){return answer;}		
};


struct border
{
	static const int radius = 1;
	uint operator() (const Matrix<uint> &in) const{
		bool hlp = false;
		hlp |= (in(1, 0) == 0) | (in(1, 2) == 0);
		hlp |= (in(0, 1) == 0) | (in(2, 1) == 0);
		if (hlp && in(1,1) == 1){
			return 1;
		}
		else{
			return 0;
		}
	}
};

typedef unsigned long long ulonglong; 

class Object
{
private:
	uint xstart = 0, xend = 0, ystart = 0, yend = 0;
	uint objid;
	bool valid = false;
	uint id = 0;
	tuple <uint, uint> mass_center, center = make_tuple(0,0);
	object_t type = GEAR;
public:
	friend double distance_beetwen_objects(const Object& ,const Object& );
	friend bool suitable_gear(const Object&, const Object&, const Object &);
	double rmin = MyINFINITY, rmax = 0;
	static const int radius = 0;
	
	void mark(Image &in){
		uint x, y;
		tie(x, y) = center;
		in(x, y) = make_tuple(255, 255, 255);
	}
	
	uint operator()(const Matrix<uint> &in) const {
		return in(0,0) == objid;
	}
	Object(uint _objid, uint x, uint y):objid(_objid), mass_center(make_tuple(0,0)){
		xstart = xend = x;
		ystart = yend = y;
	}
	Matrix<uint> GetBinPic(const Matrix<uint> &pic)
	{
		Matrix<uint> ans = (pic.submatrix(xstart, ystart, xend - xstart, yend - ystart)).deep_copy();
		for (uint i = 0; i < ans.n_rows; i++){
			for (uint j = 0; j < ans.n_cols; j++){
				ans(i,j) = uint(pic(xstart + i, ystart + j) == objid);
			}
		}
		return ans;
	}
	
	Matrix<uint> get_teeth(const Matrix<uint> &pic)
	{
		double r, x, y;
		tie(x, y) = center;
		auto binpic = GetBinPic(pic);
		double thrs = ( rmax + rmin ) / 2;
		for (uint i = 0; i < binpic.n_rows; i++){
			for (uint j = 0; j < binpic.n_cols; j++){
				r = sqrt(((xstart + i) - x) * ((xstart + i) - x) + ((ystart + j) - y) * ((ystart + j) - y));
				if ( r  < thrs ){
					binpic(i,j) = 0;
				}
			}
		}
		return binpic;
	}
	
	double compute_radius(Matrix<uint> &pic){
		double r;
		uint x, y;
		tie (x, y) = center;
		auto tmp =(pic.submatrix(xstart, ystart, xend - xstart, yend - ystart)).unary_map(*this);
		auto borderpic = tmp.unary_map(border());
		for (uint i = 1; i < borderpic.n_rows - 1; i++){
			for (uint j = 1; j < borderpic.n_cols - 1; j++){
				if (borderpic(i, j) == 1){
					r = sqrt(((xstart + i) - x) * ((xstart + i) - x) + ((ystart + j) - y) * ((ystart + j) - y));
					if (r > rmax){
						rmax = r;
					}
					if (r < rmin){
						rmin = r;
					}
				}
			}
		}
		return rmax;
		
	}
	object_t detect_type(){
		if ((rmax - rmin) < 1.5){
			return AXIS;
		}else{
			return GEAR;
		}
	}
	tuple<uint, uint> detect_mass_center(Matrix<uint> &pic){
		ulonglong mass = 0;
		ulonglong x = 0, y = 0;
		for(uint i = xstart; i <= xend; i++){
			for (uint j = ystart; j <= yend; j++){
				 mass += (pic(i,j) == id);
				  x += i * (pic(i,j) == id);
				  y += j * (pic(i,j) == id);
			  }
		  }
		  x = x / mass;
		  y = y / mass;
		  mass_center = make_tuple(x, y);
		  return make_tuple(y, x);
	}
	tuple<uint, uint> get_mass_center(){return mass_center;}
	tuple<uint, uint, uint, uint> get_bound()
	{
		return make_tuple(xstart, xend, ystart, yend);
	}
	void set_id(uint _id){
		id = _id;
		valid = true;
	}
	bool validable(){return valid;}
	void merge(Object &obj){
		if (xstart > obj.xstart)
			xstart = obj.xstart;
		if (xend < obj.xend)
			xend = obj.xend;
		if (ystart > obj.ystart)
			ystart = obj.ystart;
		if (yend < obj.yend)
			yend = obj.yend;
	}
	void update_axis(uint x, uint y){
		if (x < xstart)
			xstart = x;
		if (x > xend)
			xend = x;
		if (y < ystart)
			ystart = y;
		if (y > yend)
			yend = y;
	}
	tuple <uint, uint> center_byDT(const Matrix<uint> &m){
		uint x, y;
		uint max = 0;
		for (uint i = 0; i < m.n_rows; i++){
			for (uint j = 0; j < m.n_cols; j++){
				if (m(i,j) > max){
					x = i;
					y = j;
					max = m(i,j);
				}
			}
		}
		center = make_tuple(xstart + x, ystart + y);
		return make_tuple(ystart + y, xstart + x);
	}
	void insert(Image &in, const Image &from, Matrix<uint> grind, const Object &where)
	{
		uint x, y, x1, y1;
		tie(x, y) = where.center;
		tie(x1, y1) = center;
		uint tmpx = x - x1, tmpy = y - y1;
		for (uint i = 0; i < grind.n_rows; i++){
			for (uint j = 0; j < grind.n_cols; j++){
				if (grind(i,j) == objid){
					in(tmpx + i, tmpy + j) = from(i, j);
				}
			}
		}
	}
};

class Diversed_object
{
private: 
	uint obj_count = 0;
	uint next_obj_number = 2;
	Matrix<uint> pic;
	vector<uint> equal;
	vector<Object> obj;
	void connect(uint n, uint m);
	uint add_complies(uint n, uint m, uint i, uint j);
	uint add_obj(uint x, uint y);
    uint perfom_obj(uint n, uint x, uint y);
	void perfom();
	uint next_getting_object = 0;
public:
	vector<Object> &get_objects(){
		return obj;
	}
	uint get_obj_count() { return obj_count;}
	Diversed_object(Matrix<uint> &binpic);
	tuple<uint, uint, uint, uint> get_obj_bound(){
		if (next_getting_object == obj_count) return make_tuple(0,0,0,0);
		return obj[next_getting_object++].get_bound();
	}
	void Diver();
	Matrix<uint> &get_Matrix(){return pic;}

};


Diversed_object::Diversed_object(Matrix<uint> &binpic):pic(binpic.deep_copy()), equal(2, 0), obj()
{
	Diver();
	perfom();
}

void Diversed_object::perfom(){
	for (uint i = 2; i < next_obj_number; i++){
		if (equal[i] != i){
			equal[i] = equal[equal[i]]; // каждое элемент массива будет содеражть номер истинного объекта пиксела
			obj[equal[i] - 2].merge(obj[i - 2]);
		}
		else
		{
			uint xstart, xend, ystart, yend;
			tie(xstart, xend, ystart, yend) = obj[i - 2].get_bound();
			if (xend - xstart + yend - ystart <= 1) continue;
			obj_count++;
			obj[i - 2].set_id(i);
		}
	}
	
	for (auto i = obj.begin(); i != obj.end();){
		if (! (*i).validable()){
				obj.erase(i);
		}
		else
		{
				i++;
		}
	}
	for (uint i = 0; i < pic.n_rows; ++i)
		for (uint j = 0; j < pic.n_cols; ++j)
			pic(i,j) = equal[pic(i,j)];
			
}

void Diversed_object::connect(uint n, uint m)
{
	if( n != equal[n]){
		connect(equal[n], m);
	}
	if ( m != equal[m]){
		connect(equal[m], n);
	}
	if ( equal[m] == equal[n] )
		return;
	if (m > n){
		equal[m] = equal[n];
	}
	else{
		equal[n] = equal[m];
	}
}

uint Diversed_object::add_complies(uint n, uint m, uint newx, uint newy)
{ //добавим, что n и m обозначают одну и ту же область
	pic(newx, newy) = n;
	if ( n == m ){
		 return n;
	}
	connect(n,m);
	return n;
}

uint Diversed_object::add_obj(uint x, uint y)
{
	equal.push_back(next_obj_number);
	obj.push_back(Object(next_obj_number ,x, y));
	pic(x, y) = next_obj_number++;
	return next_obj_number;
}

uint Diversed_object::perfom_obj(uint obj_number, uint newx, uint newy)
{
	pic(newx,newy) = obj_number;
	obj[obj_number - 2].update_axis(newx, newy);
	return obj_number;
}

void Diversed_object::Diver()
{
	next_obj_number = 2;
	if (pic(0,0) != 0) {
		add_obj(0,0);
	}
	for (uint j = 1; j < pic.n_cols; j++)
	{ // пройдем по первой строчке и выделим объекты
		if (pic(0, j) != 0)
		{
			pic(0, j - 1) != 0 ? perfom_obj(pic(0, j - 1), 0, j) : add_obj(0, j); // продолжается предыдущий объект или начинается новый
		}
	}
	for (uint i = 1; i < pic.n_rows; i++)
	{
		if (pic(i, 0) != 0) // отдельно обработать первый столбец
		{
			(pic(i - 1, 0) != 0) ? perfom_obj(pic(i - 1, 0), i, 0) : add_obj(i, 0);
		}
		for (uint j = 1; j < pic.n_cols; j++)
		{
			if (pic(i, j) != 0)
			{
				if (pic(i - 1, j) == 0 && pic(i, j - 1) == 0)
				{
					add_obj(i ,j);
					continue;
				}
				if (pic(i - 1,j) != pic(i, j - 1) && pic(i - 1, j) != 0 && pic(i, j - 1) != 0)
				{
					add_complies(pic(i - 1, j), pic(i, j - 1), i, j);
					continue;
				}
				 // остался случай, когда один из pic() равен нулю, а другой нет
				 
				 if (pic(i, j - 1) != 0 && pic(i - 1, j) == 0){
					perfom_obj(pic(i, j - 1), i, j);
					continue;
				}
				if (pic(i - 1, j) != 0 && pic(i, j - 1) == 0){
					perfom_obj(pic(i - 1, j), i, j);
					continue;
				}
				perfom_obj(pic(i - 1, j), i , j);
			}
		}
	}
}

bool suitable_gear(const Object &q, const Object &one, const Object &where)
{
	uint a, b, c, d;
	tie(a,b) = where.center;
	tie(c, d) = one.center;
	double r = sqrt( (a - c) * (a - c) + (b - d) * (b - d));
	return (r > q.rmax + one.rmin && r > q.rmin + one.rmax && r < q.rmax + one.rmax);
}

double distance_beetwen_objects(const Object& one, const Object& two){
	uint x1, x2, y1, y2;
	tie(x1, y1) = one.center;
	tie(x2, y2) = two.center;
	return sqrt( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

uint max(uint a, uint b, uint c){
	if (a > b){
		return (a > c)? a : c;
	}
	else{
		return (b > c)?  b : c;
	}
}

template <typename ValueT>
class Histogram
{
private:
	const uint size;
	vector<ValueT> hist;
public:
	Histogram(uint _size, ValueT value);
	~Histogram ();
	ValueT sum() const;
	ValueT &operator [](uint i);
	double math_expection();
};

template <typename ValueT>
Histogram<ValueT>::Histogram(uint _size, ValueT value): size(_size), hist(_size, value){}	
		
template<typename ValueT>
Histogram<ValueT>::~Histogram(){}
		
template<typename ValueT>
ValueT Histogram<ValueT>:: sum() const
{
	ValueT ans = 0;
	for (auto i = hist.begin(); i != hist.end(); i++)
		ans += *i;
	return ans;
}

template<typename ValueT>
ValueT &Histogram<ValueT>:: operator [](uint i)
{
	return hist[i];
}

template<typename ValueT>
double Histogram<ValueT>:: math_expection()
{
	ValueT ans = 0;
	for (uint i = 0; i < size; i++)
	{
		ans += i * hist[i];
	}
	return ans;
}


class rgb2bright
{
public:
	Histogram<uint> hist = Histogram<uint>(256, 0); //create histogram with 256 elements, which have value 0
	uint operator() (const Image &im)
	{
		uint r,g,b,l;
		tie(r,g,b) = im(0,0); 
		l = max(r,g,b);
		hist[l]++;
		return l;
	}
	uint otsu_threshold()
	{
		uint allsum = hist.sum();
		Histogram<double> probably_hist(256, 0);
		for (uint i = 0; i < 256; i++)
		{
			 probably_hist[i] = double(hist[i]) / allsum;
		}
		double p = 0, exp = probably_hist.math_expection();
		double part_expection = 0;
		double m1, m2;
		double disp, dispmax = 0;
		uint thershold = 0;
		for (uint i = 0; i < 255; i++)
		{
			p += probably_hist[i];
			part_expection += i * probably_hist[i];
			m1 = part_expection / p;
			m2 = (exp - part_expection) / (1 - p);
			disp = p * (1 - p) * (m2 - m1)*(m2 - m1);
			if (disp > dispmax)
			{
				thershold = i;
				dispmax = disp;
			}
		}
		return thershold;
	};
	static const int radius = 0;
}; 

class bright2bin
{
public:
	uint threshold;
	static const int radius = 0;
	bright2bin(uint thrs): threshold(thrs){}
	uint operator ()(const Matrix<uint> &im) const{
		return im(0,0) > threshold ? 1 : 0;
	}
};

tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in, const char* picname)
{
    // Base: return array of found objects and index of the correct gear
    // Bonus: return additional parameters of gears
    auto object_array = vector<shared_ptr<IObject>>();
	rgb2bright obj_for_hist;
	auto Imag = in.unary_map(obj_for_hist);
	auto ans = Imag.unary_map(bright2bin(obj_for_hist.otsu_threshold()));
	Diversed_object pic(ans);
	auto obj = pic.get_objects();
	auto iter = obj.begin();
	for (auto i = obj.begin(); i != obj.end(); i++)
	{
		uint a, b, c, d;
		shared_ptr<IObject> tmp;
		bool is_broken;
		DistanceTransform q((*i).GetBinPic(pic.get_Matrix()));
		q.distance_transform();
		tie(c, d) = (*i).center_byDT(q.get_answer());
		(*i).compute_radius(pic.get_Matrix());
						
		if ( (*i).detect_type() == GEAR)
		{
				tie(a, b) = (*i).detect_mass_center(pic.get_Matrix());
				auto teeth = (*i).get_teeth(pic.get_Matrix());
				Diversed_object hlp(teeth);
				int obj_count = hlp.get_obj_count();
				is_broken =  ((a - c) * (a - c) + (b - d) * (b - d)) > 8;
				tmp.reset(new Gear(make_tuple(c, d), (*i).rmin, (*i).rmax, is_broken, obj_count));
				object_array.push_back(tmp);
				if (is_broken){
					iter = i;
				}
		} else
		{
				tmp.reset( new Axis(make_tuple(c,d)));
				object_array.push_back(tmp);
				iter = i;
		}
	}
	
	double dist1 = MyINFINITY, dist2 = MyINFINITY;
	
	auto one = obj.begin(), two = obj.begin();
	
	
	for (auto i = obj.begin(); i != obj.end(); i++){
		if (i == iter)
			continue;
		double tmp = distance_beetwen_objects(*i, *iter);
		if (tmp < dist2){
			if (tmp < dist1){
				dist2 = dist1;
				two = one;
				dist1 = tmp;
				one = i;
			}else{
				dist2 = tmp;
				two = i;
			}
		}
	}
	
		
	Image answer_image = in.deep_copy();
    int result_idx = 41;
	for (int i = 1; i <= 3; i++){
		char *gearname = new char[strlen(picname) + 3];
		uint j;
		sprintf(gearname, "%s", picname);
		for (j = strlen(picname) - 1; j != 0 && picname[j] != '.'; j--);
		if (j == 0) throw "mistake in name of image";
		gearname[j] = '\0';
		sprintf(gearname, "%s_%d.bmp", gearname, i);
		Image detail = load_image(gearname);
		rgb2bright otherobjforhist;
		auto tmpobj = detail.unary_map(otherobjforhist);
		auto binaryimg = tmpobj.unary_map(bright2bin(otherobjforhist.otsu_threshold()));
		Diversed_object gear(binaryimg);
		if (gear.get_obj_count() != 1) throw "invalid gear image";
		auto gear_obj = gear.get_objects();
		DistanceTransform qw(gear_obj[0].GetBinPic(gear.get_Matrix()));
		qw.distance_transform();
		gear_obj[0].center_byDT(qw.get_answer());
		gear_obj[0].compute_radius(gear.get_Matrix());
		if (suitable_gear(gear_obj[0], *one, *iter) && suitable_gear(gear_obj[0], *two, *iter)){
			gear_obj[0].insert(answer_image, detail, gear.get_Matrix(), *iter);
			result_idx = i;
			break;
		}
		delete []gearname;
	}


    return make_tuple(result_idx, object_array, answer_image);

}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        ofstream fout(argv[3]);

        vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;
        tie(result_idx, object_array, dst_image) = repair_mechanism(src_image,  argv[1]);
        save_image(dst_image, argv[2]);

        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
