
# IniFileParser Class

[[_TOC_]]


## Introduction

The `IniFileParser` class provides a simple and efficient way to read and parse INI configuration files in C++. It supports reading values from different sections and keys within the INI file and also offers utility methods to retrieve integer and double values. This class includes the following features:
- Parse sections and key-value pairs from an INI file.
- Retrieve values as `string`, `int`, and `double`.
- Get sections with a specific prefix.
- Handles various common INI file formatting styles (comments, empty lines, etc.).


### Example INI File Format

```
[Settings]
LogLevel = debug
MaxRetries = 5

[Database]
Host = localhost
Port = 5432
MaxConnections = 100

[Network]
Timeout = 30.5
RetryInterval = 5.0
```

### Constructor

The constructor `IniFileParser(const std::string& filenameAux)` takes a file name as input and loads the .ini file, storing its contents in a structured format. Here’s how it works:
1. File Opening: The constructor attempts to open the specified file. If the file cannot be opened, it displays an error message and exits.
2. Section and Key-Value Pair Extraction:
    - Each line of the file is read and trimmed of whitespace.
    - Lines beginning with `[` are recognized as section headers, with the section name extracted and added to an internal `sections` map. A separate `sectionOrder` vector also records section names to maintain the order as they appear in the file.
    - Lines containing `=` are processed as key-value pairs. The part before `=` is the key, and the part after is the value, both of which are trimmed and added under the current section.
3. Comments and Malformed Lines: Lines starting with `;` or `#` are treated as comments and skipped. Lines missing an expected `=` or `]` are flagged as malformed and ignored.

This parser stores the .ini file data in a map of maps, sections, for easy access by section and key. Additionally, `sectionOrder` ensures sections are retrieved in the file’s original order, which can be useful for preserving structure in applications where order matters.


### Methods

#### `std::string getValue(const std::string& section, const std::string& key) const`

This method retrieves the value of a specific key from a given section of the INI file. It returns an empty string if the section or key does not exist.

**Parameters:**
- `section`: The section in the INI file.
- `key`: The key within the section.

**Example:**
```cpp
std::string value = parser.getValue("Settings", "LogLevel");
```

#### `std::vector<std::map<std::string, std::string>> getSectionsData(const std::string& sectionPrefix) const`

This method retrieves all sections that start with a specific prefix. It returns a `vector` containing the key-value pairs from all matching sections.

**Parameters:**
- `sectionPrefix`: The prefix to match the section names.

**Example:**
```cpp
std::vector<std::map<std::string, std::string>> sections = parser.getSectionsData("Database");
```

#### `std::vector<std::pair<std::string, std::map<std::string, std::string>>> getSections(const std::string& sectionPrefix) const`

This method retrieves all sections that start with a specific prefix. It returns a `vector` containing a pair of values, a string with the name of the section and a second, with key-value pairs from all matching sections.

**Parameters:**
- `sectionPrefix`: The prefix to match the section names.

**Example:**
```cpp
std::vector<std::pair<std::string, std::map<std::string, std::string>>> thermodynamicData = parser.getSections("Data");
```

#### `int getInt(const std::string& section, const std::string& key) const`

This method retrieves an integer value from a specific section and key of the INI file. It throws an exception if the value cannot be converted to an integer.

**Parameters:**
- `section`: The section in the INI file.
- `key`: The key within the section.

**Example:**
```cpp
int port = parser.getInt("Database", "Port");
```

#### `int getInt(const std::string& section, const std::string& key, int defaultValue) const`

This overloaded method retrieves an integer value associated with a specific key in a given section of the INI file. If the key is missing or its value is empty, the method returns the specified `defaultValue`. If the key exists but its value cannot be converted to an integer, the method may throw an exception (e.g., std::invalid_argument).

**Parameters:**
- `section`: The section in the INI file.
- `key`: The key within the section.
- `defaultValue`: The value to return if the key is missing or its value is empty.

**Example:**
```cpp
int port = parser.getInt("Database", "Port", 20);
```

#### `int getInt(const std::map<std::string, std::string>& section, const std::string& key) const`

This overloaded method retrieves an integer value from a map representing a section, instead of looking it up from the entire INI file.

**Parameters:**
- `section`: A map representing a section.
- `key`: The key within the map.

**Example:**
```cpp
std::map<std::string, std::string> dbSettings = parser.getSectionsData("Database")[0];
int maxConnections = parser.getInt(dbSettings, "MaxConnections");
```

#### `double getDouble(const std::string& section, const std::string& key) const`

This method retrieves a double value from a specific section and key. It throws an exception if the value cannot be converted to a double.

**Parameters:**
- `section`: The section in the INI file.
- `key`: The key within the section.

**Example:**
```cpp
double timeout = parser.getDouble("Network", "Timeout");
```

#### `double getDouble(const std::string& section, const std::string& key, double defaultValue) const `

This overloaded method retrieves a double value associated with a specific key in a given section of the INI file. If the key is missing or its value is empty, the method returns the specified `defaultValue`. If the key exists but its value cannot be converted to a double, the method may throw an exception (e.g., std::invalid_argument).

**Parameters:**
- `section`: The section in the INI file.
- `key`: The key within the section.
- `defaultValue`: The value to return if the key is missing or its value is empty.

**Example:**
```cpp
double timeout = parser.getDouble("Network", "Timeout", 50.5);
```

#### `double getDouble(const std::map<std::string, std::string>& section, const std::string& key) const`

This overloaded method retrieves a double value from a map representing a section.

**Parameters:**
- `section`: A map representing a section.
- `key`: The key within the map.

**Example:**
```cpp
std::map<std::string, std::string> netSettings = parser.getSectionsData("Network")[0];
double retryInterval = parser.getDouble(netSettings, "RetryInterval");
```

#### `static std::string trim(const std::string& str)`

This helper method removes leading and trailing whitespace characters from a string.

**Parameters:**
- `str`: The string to be trimmed.


## Usage

### Include the Header

```cpp
#include "ini_file_parser/IniFileParser.h"
```

### Create a `IniFileParser` Object

You can instantiate the `IniFileParser` by providing the path to the INI file you want to parse:

```cpp
IniFileParser parser("config.ini");
```
