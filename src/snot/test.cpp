#include <iostream>


#include "param_storage.h"
#include "snot.h"

using namespace std;
using namespace gsparams;

int main(){


    string testin="{\"tfs\":{\"A\":[1.1,1.2,1.3],\"B\":[2.1,2.2,2.3],\"C\":[2.1,2.2,2.3,2.1,2.2,2.3]}}";

    DictList my_dictlist = parse_dictlist(testin);


        std::vector<double> thevals(0);


        for(int j = 0;j<5;j++){
        my_dictlist.traverse(&thevals);
        for(int i = 0;i<thevals.size();i++){
            cout << thevals[i] << " ";
        }
        cout << endl;
        }
        thevals.clear();


        cout << my_dictlist.str() << endl;
    return 0;
}
