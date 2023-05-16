/*
 * LND.cpp
 *
 *  Created on: Apr 26, 2023
 *      Author: lyx
 */
#include <iostream>
#include <string>
#include <array>
#include <vector>

using namespace std;

vector<int> load_image(string imagefile)
{
	cout<<imagefile;
    vector<int> imagedata(100,0);
    return imagedata;
}

class FM
{
private:
	string imgfile;
	vector<int> imgdata;
public:
	FM (string imagefile)
	{
		imgfile = imagefile;
		imgdata = load_image(imagefile);
	}
};

int main()
{
	string f;
	cin>>f;
	cout<<"LND project is for local neurites detector\n";
	cout<<"along with soma detection\n";
	cout<<"\nNOW it is doing "<<f<<"\n\n";

	string imgfile;
	cin>>imgfile;
	vector<int> imgdata;
	imgdata = load_image(imgfile);
	return 0;
}
