#ifndef GS_PARSING_H
#define GS_PARSING_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>


class SimpleStringTokenizer{
	public:
		SimpleStringTokenizer(){};

		int tokenize(std::string in){
			//https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
			thestring = in;
			tokens.clear();
			line_ss.clear();
			line_ss.str(in);
			copy(std::istream_iterator<std::string>(line_ss),
				std::istream_iterator<std::string>(),
				std::back_inserter(tokens));
			return tokens.size();
		}

		std::string get_string() const{
			return thestring;
		}

		const std::vector< std::string >& get_tokens() const{
			return tokens;
		};

		const std::string& operator[]( std::vector<std::string>::size_type pos ) const{
			return tokens[pos];
		}

		int size() const{
			return tokens.size();
		}
	private:
		std::string thestring;
		std::istringstream line_ss;
		std::vector< std::string > tokens;
};

#endif
