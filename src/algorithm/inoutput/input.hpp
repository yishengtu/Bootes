#ifndef INPUT_HPP_
#define INPUT_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <typeinfo>


using namespace std;

class input_file{
    private:
        string fn_;
    public:
        map<string, string> inputdict;
    input_file(string fn){
        fn_ = fn;
        ifstream infile(fn);
        for (string line; getline(infile, line); ){
            int comment_ind = line.find('#');
            string pure_info = line.substr(0, comment_ind);
            int equsign = line.find('=');
            if (equsign == -1){
                continue;   // If "=" is not found, then skip this line
            }
            string var_name = pure_info.substr(0, equsign);
            string var_name_clean = clean(var_name);
            string var_value = pure_info.substr(equsign + 1, -1);
            string var_value_clean = clean(var_value);

            inputdict.insert({var_name_clean, var_value_clean});
        }
    }

    int getInt(string name){
        try{
            return stoi(inputdict[name]);
        }
        catch (std::invalid_argument e){
            cout << "\"" << name << "\"" << " does not exist in input file" << endl << flush;
            throw e;
        }
    }

    float getFloat(string name){
        try{
            return stof(inputdict[name]);
        }
        catch (std::invalid_argument e){
            cout << "\"" << name << "\"" << " does not exist in input file" << endl << flush;
            throw e;
        }
    }

    double getDouble(string name){
        try{
            return stod(inputdict[name]);
        }
        catch (std::invalid_argument e){
            cout << "\"" << name << "\"" << " does not exist in input file" << endl << flush;
            throw e;
        }
    }

    string getString(string name){
        try{
            return inputdict[name];
        }
        catch (std::invalid_argument e){
            cout << "\"" << name << "\"" << " does not exist in input file" << endl << flush;
            throw e;
        }
    }

    string clean(string var){
        string var_out = var;
        /* remove all space before value */
        int first_spaceIND = var_out.find(' ');
        while (first_spaceIND == 0){
            var_out = var_out.substr(1, -1);
            first_spaceIND = var_out.find(' ');
        }
        /*remove all space after value */
        int last_spaceIND = var_out.find_last_of(' ');
        while (last_spaceIND == var_out.length() - 1){
            var_out = var_out.substr(0, last_spaceIND);
            last_spaceIND = var_out.find_last_of(' ');
        }
        return var_out;
    }
};

#endif
