#ifndef GS_PARAM_STORAGE
#define GS_PARAM_STORAGE

/*
MIT License

Copyright (c) 2017 Bryan Lunt <bjlunt2@illinois.edu> <bryan.j.lunt@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <stdexcept>
#include <vector>
#include <map>
#include <string>
#include <stack>
#include <sstream>

#include <iostream>

namespace gsparams {

#define dictlist_primitive_t double
#define dictlist_key_t std::string

typedef enum { undecided, primitive, dict, list} DictListType;

class DictList {
    public:
    //private:
    //public:
        DictListType my_type;
        dictlist_primitive_t my_value;

        typedef std::map<std::string,int> dictlist_map_t;
        dictlist_map_t map_storage;
        std::vector<std::string> map_key_storage; //when we want to output, we need to be able to iterate over keys too.
        std::vector<DictList> list_storage; //only primitive values lack list storage.

        inline void undecided_to_dict_else_error(){
            if(this->my_type == dict){return;
            }if(this->my_type == undecided){
                list_storage.clear();
                //list_storage = new std::vector<DictList>(0);
                map_storage.clear();
                //map_storage = new dictlist_map_t();
                map_key_storage.clear();
                //map_key_storage = new std::vector<std::string>(0);
                this->my_type = dict;
            }else{
                throw std::runtime_error("Already upgraded");
            }
        }

        inline void undecided_to_list_else_error(){
            if(this->my_type == list){return;
            }if(this->my_type == undecided){
                list_storage.clear();
                //list_storage = new std::vector<DictList>(0);
                map_storage.clear();
                //map_storage = NULL;
                map_key_storage.clear();
                //map_key_storage = NULL;
                this->my_type = list;
            }else{
                throw std::runtime_error("Already upgraded");
            }
        }

    //protected:
        inline void traverse_internal(std::vector< dictlist_primitive_t >* target) const {
            //public function has already cleared and setup the beginnigs of the target vector.
            if(undecided == this->my_type){
                return;
            }

            if(primitive == this->my_type){
                target->push_back(this->v());
                return;
            }

            int n = this->list_storage.size();
            for(int i = 0; i < n; i++){
                this->list_storage.at(i).traverse_internal(target);
                //shoudl be able to overwrite some operator with an iterator and
            }
        }

        inline int populate_internal(const std::vector< dictlist_primitive_t >& target,int starting_at) {
            //public function has already cleared and setup the beginnigs of the target vector.
            if(undecided == this->my_type){
                return 0;
            }

            if(primitive == this->my_type){
                this->my_value = target[starting_at];
                return 1;
            }

            //some type that stores things

            int num_consumed = 0;
            int n = this->list_storage.size();
            for(int i = 0; i < n; i++){
                num_consumed += this->list_storage.at(i).populate_internal(target,starting_at+num_consumed);
            }
            return num_consumed;
        }

        inline void clear_helper(){
            this->list_storage.clear();
            this->map_key_storage.clear();
            this->map_storage.clear();
            this->my_value = 0.0;
            this->my_type = undecided;
        }

        inline void copy_helper(const DictList &other){

            this->my_type = other.my_type;
            this->my_value = other.my_value;

            this->list_storage = other.list_storage;
            this->map_key_storage = other.map_key_storage;
            this->map_storage = dictlist_map_t(other.map_storage);

        }
    public:
        inline DictList() : my_type(undecided), map_storage(), map_key_storage(), list_storage(), my_value(0.0) {
            my_type = undecided;
            //map_storage = NULL;
            //list_storage = NULL;
            //map_key_storage = NULL;
            my_value = 0.0;//TODO: can we make this nan or some other defensive value?
        }
        inline DictList(dictlist_primitive_t in) : my_type(undecided), map_storage(), map_key_storage(), list_storage(), my_value(0.0)  {
            my_type = primitive;
            //map_storage = NULL;
            //list_storage = NULL;
            //map_key_storage = NULL;
            my_value = in;
        }
        //copy constructor


        inline DictList(const DictList &other) : my_type(undecided), map_storage(), map_key_storage(), list_storage(), my_value(0.0)  {
            this->copy_helper(other);
        }

        inline DictList& operator=(const DictList &other){
            this->clear_helper();

            this->copy_helper(other);
            return *this;
        }

        inline DictList& operator=(const dictlist_primitive_t &other){
            if(undecided != this->my_type && primitive != this->my_type){
                throw std::runtime_error("Attempt to assign primitive value to nonprimitive, non undecided DictList.");
            }
            this->clear_helper();
            this->my_type = primitive;
            this->my_value = other;

            return *this;
        }

        inline ~DictList(){
        }

        class iterator;

        inline dictlist_primitive_t v() const {
            if(this->my_type != primitive){throw std::runtime_error("Asked value of non primitive");}
            return this->my_value;
        }

        inline void push_back(const DictList& in) {
            switch(this->my_type){
                case primitive:
                    throw std::runtime_error("Cannot append to a primitive");
                    break;
                case undecided:
                    this->my_type = list;
                    //this->list_storage = new std::vector<DictList>(0);
                    this->list_storage.clear();
                case list:
                    this->list_storage.push_back(in);
                    break;
                case dict:
                    throw std::runtime_error("Cannot append to a dict");
                    break;
                default:
                    throw std::runtime_error("somethind strange.");
            }
        }

        inline void push_back(dictlist_primitive_t in){
            #ifdef GS_PARAM_STORAGE_DEBUG
            std::cerr << "append from value" << std::endl;
            #endif
            DictList indictlist(in);
            this->push_back(indictlist);
        }

        inline DictList& at(const int location) {
            switch(this->my_type){

                case list:
                case dict:
                    //in both of these cases, we access as a list.
                    return this->list_storage.at(location);//vector::at() is in C++98, keep
                    break;
                case primitive:
                    throw std::runtime_error("Can't subscript a primitive value.");
                    break;
                case undecided:
                default:
                    throw std::runtime_error("Cant get values from undecided dictlist.");
            }
            throw std::runtime_error("Cant get values from undecided dictlist.");
            //compiler warning
        }


        inline DictList& operator[](const int location) {
            return this->at(location);
        }


        inline void set(dictlist_key_t key, DictList& in){
            undecided_to_dict_else_error();
            if(this->my_type != dict){throw std::runtime_error("Not a dictionary.");}
            //What if the value already exists?
            dictlist_map_t::iterator key_to_int = this->map_storage.find(key);
            if(this->map_storage.end() == key_to_int){
                //key does not exist
                this->list_storage.push_back(in);
                this->map_key_storage.push_back(key);
                this->map_storage[key] = (int)(this->list_storage.size()-1);
            }else{
                //The key already exists, so we need to replace it.
                this->list_storage[key_to_int->second] = in;
            }
        }

        inline void set(dictlist_key_t key, dictlist_primitive_t in){
            DictList indictlist(in);
            this->set(key,indictlist);
        }

        inline DictList& at(dictlist_key_t key) {
            if(dict != this->my_type){
                throw std::runtime_error("Can't us this as a dictionary.");
            }
            dictlist_map_t::iterator key_to_int = this->map_storage.find(key);
            if(this->map_storage.end() == key_to_int){
                throw std::out_of_range("Asked for a key that does not exist.");
            }
            return this->list_storage.at(key_to_int->second);//vector::at() in C++98
        }

        inline DictList& operator[](dictlist_key_t key) {
            undecided_to_dict_else_error();
            dictlist_map_t::iterator key_to_int = this->map_storage.find(key);
            if(this->map_storage.end() == key_to_int){
                //
                DictList tmp;
                this->set(key,tmp);//should insert into the key list if necessary

                //Try again using the same search code.
                key_to_int = this->map_storage.find(key);
                if(this->map_storage.end() == key_to_int){
                    throw std::out_of_range("this function should have created a new one...");
                }
            }

            return this->list_storage.at(key_to_int->second);
        }

        inline DictList& operator[](const char* key){
            return this->operator[](dictlist_key_t(key));
        }

        inline int size() const {
            switch(this->my_type){
                case undecided:
                    return -2;
                    break;
                case primitive:
                    return 1;
                    break;
                case list:
                case dict:
                    return this->list_storage.size();
                default:
                    break;
            }
        }

        inline void traverse(std::vector< dictlist_primitive_t >* target) const {
            if(undecided == this->my_type){
                throw std::runtime_error("Can't traverse undecided");
            }

            target->clear();

            this->traverse_internal(target);
        }

        inline void populate(const std::vector< dictlist_primitive_t >& source) {
            int num_consumed = populate_internal(source,0);
            //do I want to check that the number consumed is right?
        }

        inline bool populate_or_revert(const std::vector< dictlist_primitive_t >& source) {
            std::vector< dictlist_primitive_t > tmp(0);
            this->traverse(&tmp);

            try{
                this->populate(source);
                return true;
            }catch(std::range_error e){
                this->populate(tmp);
                return false;
            }

            return false;
        }


class iterator : public std::forward_iterator_tag {
    protected:
        /*This iterator will need some kind of stack for state storage.
        The things it iterates over each provide iterators, so maybe the natural thing is to have a stack of iterators?
        It might not be storage efficient, but it's programmer time efficient.

        Duh. No, we don't need an explicit stack. The iterator knows which DictList it iterates, and has a handle to a current sub-iterator.
        The handles to handles to handles of sub-iterators _are_ the stack.

        Blah. That's just reasoning to avoid learning more STL, it will create a lot of memory derefrences.
        Yes, this library won't be used in really slow parts of the code, but still...
        */
        DictList* my_dictlist;
        std::stack< std::pair<DictList*, int> > my_stack;
        //std::vector<std::vector<DictList>::iterator> position_stack;

    public:
        inline iterator(DictList *initial_dictlist, int position) : my_dictlist(initial_dictlist), my_stack(){
            my_stack.push(std::make_pair(initial_dictlist,position));
        }

    public:

        inline iterator(const DictList::iterator &other) : my_dictlist(other.my_dictlist), my_stack(other.my_stack) {
            //pass
        }

        inline DictList::iterator& operator=(const DictList::iterator &other){
            my_dictlist = other.my_dictlist;
            my_stack = other.my_stack;
            return *this;
        }

        inline bool operator==(const DictList::iterator &other) const{
            return (my_dictlist == other.my_dictlist && my_stack == other.my_stack);
        }

        inline bool operator!=(const DictList::iterator &other) const{
            return !(my_dictlist == other.my_dictlist && my_stack == other.my_stack);
        }

        inline DictList::iterator& operator++(){
            //TODO: The meat of the traversal algorithm
            if(my_stack.size() < 1){ //off the end.
                //std::cerr << "off end." << std::endl;
                return *this;
            }

            bool loop_continues = false;

            while(my_stack.size() > 1 && loop_continues){
                my_stack.top().second++;
                if(my_stack.top().second >= my_stack.top().first->size()){
                    //primitives will never be on the stack, so this is ok.
                    //The stack is only those containers that we have descended into.
                    //So, I guess that implies that iterators for undecided and primitive types need to be initialized to .end();
                    my_stack.pop();
                    loop_continues = true;
                    continue;
                }

                //We're looking at the next subobject of the thing on top of the stack.
                //If it is another compound object, we need to descend into it( put it on top of the stack. )
                DictList *current_pointed_element = &(my_stack.top().first->at(my_stack.top().second));
                DictListType check_type = current_pointed_element->my_type;

                //some compound type, gets pushed to the stack, equivalent to recursion
                if(undecided != check_type && primitive != check_type){
                    my_stack.push(std::make_pair(current_pointed_element,-1));
                    continue;//need to descend into that.
                }

                if(undecided == check_type){
                    //Yes, iterators now iterate over undecided.
                    loop_continues = false;
                    break;
                }
                if(primitive == check_type){
                    loop_continues = false;
                    break;
                }

                throw std::runtime_error("Never make it here.");
            }
            /*Otherwise, we are pointing to a primitive type (or undecided, behaviour not defined)
            */


            return *this;
        }

        inline DictList::iterator operator++(int unused){
            DictList::iterator tmp(*this);
            this->operator++();
            return tmp;
        }

        inline DictList& operator*(){

            //Actually, the increment should handle this.
            DictList& current_pointed_element = my_stack.top().first->at(my_stack.top().second);

            if(primitive == current_pointed_element.my_type){
                return current_pointed_element;//TODO: need to check that this actually returns a reference, rather than creating a new one with copy construction.
            }

            //DEBUG
            //std::cerr <<


            throw std::runtime_error("Invalid iterator state.");
        }

};//End of declaration of iterator

    inline iterator begin(){
        iterator ret_iter(this,-1);
        ++ret_iter;
        return ret_iter;
    }

    inline iterator end(){
        iterator ret_iter(this,this->size());
        ++ret_iter;
        return ret_iter;
    }

    /*
    Unfortunately, it looks like we can't have typecasting and have nice subscripting at the same time.

    We can, by virtue of altering the subscript operator, which somehow makes that take presedence over this cast.
    Thanks stack overflow! : https://stackoverflow.com/questions/15850840/ambiguous-overload-for-operator-if-conversion-operator-to-int-exist

    The overloaded subscript is next to the other subscript.
    */

    inline operator dictlist_primitive_t() const {
        if(primitive != this->my_type){ throw std::runtime_error("Cannot cast non-primitive.");}
        return this->v();
    }

    friend std::ostream& operator<<(std::ostream& os, const DictList& obj);

    inline std::string str(){
        std::stringstream ss;
        ss.clear();
        ss << *this;
        std::string foobar(ss.str());
        return foobar;
    }

};//END OF DICTLIST

    inline std::ostream& operator<<(std::ostream& os, const DictList& obj)
    {
        // write obj to stream
        if(primitive == obj.my_type){
            os << obj.v();
            return os;
        }

        if(dict == obj.my_type){
            os << "{";
                int n_members = obj.size();
                if(n_members > 0){
                    os << '"' << obj.map_key_storage.at(0) << '"' << ":" << obj.list_storage.at(0);

                    for(int i = 1;i<n_members;i++){
                        os << "," << '"' << obj.map_key_storage.at(i) << '"' << ":" << obj.list_storage.at(i);
                    }

                }
            os << "}";
            return os;
        }else if(list == obj.my_type){
            os << "[";
            int n_members = obj.size();
            if(n_members > 0){
                os << obj.list_storage.at(0);

                for(int i = 1;i<n_members;i++){
                    os << "," << obj.list_storage.at(i);
                }

            }

            os << "]";
        }else if(primitive == obj.my_type){
            os << obj.v();
        }


        return os;
    }

}

#endif
