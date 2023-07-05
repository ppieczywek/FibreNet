#pragma once
#include <iostream>
#include <fstream>
#include <vector>
class TextParser {
public:
	static bool ReadFileContent(std::string inputFile, std::vector<std::string>& fileContents) {
		std::ifstream	file(inputFile);
		if (file.is_open()) {
			fileContents.clear();
			std::string str;
			while (std::getline(file, str)) fileContents.push_back(str);
			file.close();
			return true;
		}
		else {
			return false;
		}
	}

	static int ScanText(std::vector<std::string>& text, std::string input, std::vector<std::string>& tokens, int startLine) {
		std::string	 str;
		tokens.clear();
		if (startLine >= 0) {
			if (!text.empty()) {
				for (int ii = startLine; ii < static_cast<int>(text.size()); ii++) {
					str = text[ii];
					std::istringstream iss(str);
					std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
					if (!tokens.empty()) {
						if (tokens[0] == input)	return ii + 1;
					}
					tokens.clear();
				}
			}
		}
		return -1;
	}
};