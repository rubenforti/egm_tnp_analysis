#include <string>
#include <regex>
#include <boost/algorithm/string/replace.hpp>

void createjobs() {
	std::ofstream script("script.sh");
	for (unsigned int i=1; i!=18; i++) {
		std::string filename("runFits_");
		filename+=std::to_string(i)+std::string(".py");
		std::ofstream pythonfile;
		pythonfile.open(filename.c_str());
		std::ifstream infile;
		infile.open("runFits.py.tmpl");
		std::string line;
		while (std::getline(infile,line)) {
			boost::replace_all(line, "{***}", std::to_string(i).c_str());
			pythonfile<<line;
			pythonfile<<"\n";
		}
		pythonfile.close();
		infile.close();
		script<<"python ";
		script<<filename.c_str();
		script<<" -o bin";
		script<<i;
		script<<"\n";
	}
	script.close();
}
