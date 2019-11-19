/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef CXXARGS_CXXARGS_HPP
#define CXXARGS_CXXARGS_HPP

#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <utility>
#include <exception>
#include <memory>

#define CXXARGS_VERSION_MAJOR 1
#define CXXARGS_VERSION_MINOR 1
#define CXXARGS_VERSION_PATCH 1

namespace cxxargs {
  namespace exceptions {
    struct cxxargs_exception : public std::exception {
      std::string msg;
      virtual ~cxxargs_exception() override;
      cxxargs_exception(cxxargs_exception&&) = default;
      cxxargs_exception(std::string m) : msg(m) {}
      const char* what() const noexcept override { return msg.c_str(); }
    };
    struct argument_uninitialized_exception : public cxxargs_exception {
      argument_uninitialized_exception(std::string name)
	: cxxargs_exception("Argument " + name + " was not given and has no default value.") {}
      argument_uninitialized_exception(char name) : argument_uninitialized_exception(std::string(1, name)) {}
    };
    struct argument_not_defined_exception : public cxxargs_exception {
      argument_not_defined_exception(std::string name)
	: cxxargs_exception("Argument " + name + " is not defined.") {}
      argument_not_defined_exception(char name) : argument_not_defined_exception(std::string(1, name)) {}
    };
    cxxargs_exception::~cxxargs_exception() = default;
  }
  template <typename T> std::istream& operator>> (std::istream &in, std::vector<T> &t) {
    std::string str;
    while (getline(in, str, ',')) {
      t.emplace_back(*(new T()));
      std::stringstream sstream(str);
      sstream >> t.back();
    }
    return in;
  }

  class Argument {
   private:
    std::string help_text;

   public:
    Argument();
    Argument(char short_name, std::string h_text) : help_text('-' + std::string(1, short_name) + '\t' + h_text) {}
    Argument(std::string long_name, std::string h_text) : help_text(long_name + '\t' + h_text) {}
    Argument(char short_name, std::string long_name, std::string h_text) : help_text("-" + std::string(1, short_name) + ' ' + long_name + "\t" + h_text) {}
    virtual ~Argument();

    virtual void parse_argument(std::stringstream &str) =0;
    virtual void parse_argument(std::vector<std::string>::const_iterator iter) =0;
    template <class T, class U> void set_val(U& in_val);

    virtual const bool& is_initialized() const =0;
    virtual const bool& is_required() const =0;
    template <class T> const T& get_val() const;
    const std::string& get_help() const { return this->help_text; }
    virtual void set_not_required() =0;
  };
  Argument::~Argument() = default;

  template <typename T>
  class ArgumentVal : public Argument {
   private:
    T val;
    bool value_initialized = false;
    bool required = true;

   public:
    using Argument::Argument;
    ArgumentVal(char short_name, std::string help_text, T in_val)
      : Argument(short_name, help_text) {
      this->set_val(in_val);
    }
    ArgumentVal(std::string long_name, std::string help_text, T in_val)
      : Argument(long_name, help_text) {
      this->set_val(in_val);
    }
    ArgumentVal(char short_name, std::string long_name, std::string help_text, T in_val)
      : Argument(short_name, long_name, help_text) {
      this->set_val(in_val);
    }
    ~ArgumentVal() override = default;

    void parse_argument(std::stringstream &str) override {
      T in_val;
      str >> in_val;
      this->set_val(in_val);
    }
    void parse_argument(std::vector<std::string>::const_iterator iter) override {
      ++iter;
      std::stringstream str(*iter);
      parse_argument(str);
    }
    void set_val(T& in_val) { this->value_initialized = true; this->val = in_val; }

    const bool& is_initialized() const override { return this->value_initialized; }
    const bool& is_required() const override { return this->required; }
    const T& get_val() const { return this->val; }
    void set_not_required() override { this->required = false; }
  };
  template<> void ArgumentVal<bool>::parse_argument(std::vector<std::string>::const_iterator) {
    bool in_val = (this->is_initialized() ? !this->get_val() : true);
    this->set_val(in_val);
  }
  template<class T, class U> void Argument::set_val(U& in_val) {
    return dynamic_cast<ArgumentVal<T>&>(*this).set_val(in_val);
  }
  template<class T> const T& Argument::get_val() const {
    return dynamic_cast<const ArgumentVal<T>&>(*this).get_val();
  }

  class Arguments {
   private:
    std::map<std::string, std::shared_ptr<Argument>> longargs;
    std::map<char, std::shared_ptr<Argument>> shortargs;
    std::vector<std::string> positionals;
    std::string help_text;
    std::string program_name;

    template <typename T> void validate(const std::map<T, std::shared_ptr<Argument>> &args) const {
      for (auto kv : args) {
	if (!kv.second->is_initialized() && kv.second->is_required()) {
	  throw exceptions::argument_uninitialized_exception(kv.first);
	}
      }
    }
    const std::shared_ptr<Argument>& get_val(const std::string &name) const { return this->longargs.at(name); }
    const std::shared_ptr<Argument>& get_val(const char &name) const { return this->shortargs.at(name); }
    template<typename T> void set_own_val(const std::string &name, T in_val) { this->longargs.at(name)->set_val<T>(in_val); }
    template<typename T> void set_own_val(const char &name, T in_val) { this->shortargs.at(name)->set_val<T>(in_val); }

   public:
    Arguments(std::string p_name, std::string u_info)
      : help_text(u_info), program_name(p_name) {}

    template <typename T> void set_val(const char &name, T in_val) {
      this->set_own_val<T>(name, in_val);
    }
    template <typename T> void set_val(const std::string& name, T in_val) {
      this->set_own_val<T>("--" + name, in_val);
    }

    template <typename T> void add_argument(char s_name, std::string l_name, std::string h_text) {
      this->longargs.insert(std::make_pair("--" + l_name, std::shared_ptr<Argument>(new ArgumentVal<T>(s_name, "--" + l_name, h_text))));
      this->shortargs.insert(std::make_pair(s_name, this->longargs.at("--" + l_name)));
      this->help_text += '\n' + this->longargs.at("--" + l_name)->get_help();
    }
    template <typename T> void add_argument(char s_name, std::string l_name, std::string h_text, T in_val) {
      this->add_argument<T>(s_name, l_name, h_text);
      this->set_val<T>(s_name, in_val);
      this->set_val<T>(l_name, in_val);
    }

    template <typename T> void add_long_argument(std::string l_name, std::string h_text) {
      this->longargs.insert(std::make_pair("--" + l_name, std::shared_ptr<Argument>(new ArgumentVal<T>("--" + l_name, h_text))));
      this->help_text += '\n' + this->longargs.at("--" + l_name)->get_help();
    }
    template<typename T> void add_long_argument(std::string l_name, std::string h_text, T in_val) {
      this->add_long_argument<T>(l_name, h_text);
      this->set_val<T>(l_name, in_val);
    }

    template <typename T> void add_short_argument(char s_name, std::string h_text) {
      this->shortargs.insert(std::make_pair(s_name, std::shared_ptr<Argument>(new ArgumentVal<T>(s_name, h_text))));
      this->help_text += '\n' + this->shortargs.at(s_name)->get_help();
    }
    template <typename T> void add_short_argument(char s_name, std::string h_text, T in_val) {
      this->add_short_argument<T>(s_name, h_text);
      this->set_val(s_name, in_val);
    }

    void set_not_required(const char &s_name) {
      this->shortargs.at(s_name)->set_not_required();
    }
    void set_not_required(const std::string &l_name) {
      this->longargs.at("--" + l_name)->set_not_required();
    }

    void parse(int argc, char** argv) {
      std::vector<std::string> vec(argv, argv+argc);      
      for (std::vector<std::string>::const_iterator it = vec.begin() + 1; it < vec.end(); ++it) {
	if (this->longargs.find(*it) != this->longargs.end()) {
	  this->longargs.at(*it)->parse_argument(it);
	} else if (it->compare(0, 1, "-") == 0 && it->compare(1, 1, "-") != 0) {
	  for (size_t i = 1; i < it->size(); ++i) {
	    if (this->shortargs.find(it->at(i)) != this->shortargs.end()) {
	      this->shortargs.at(it->at(i))->parse_argument(it);
	    }
	  }
	} else if (it->find('=') != std::string::npos) {
	  std::stringstream arg(*it);
	  std::string name;
	  getline(arg, name, '=');
	  if (this->longargs.find(name) != this->longargs.end()) {
	    this->longargs.at(name)->parse_argument(arg);
	  }
	} else if (it->compare("--") == 0) {
	  while (++it < vec.end()) {
	    this->positionals.emplace_back(*it);
	  }
	}
      }
      this->validate(this->longargs);
      this->validate(this->shortargs);
    }

    size_t n_positionals() const { return this->positionals.size(); }

    template <typename T> const T& value(const char &name) const {
      return (*this->get_val(name)).template get_val<T>();
    }
    template <typename T> const T& value(const std::string &name) const {
      return (*this->get_val("--" + name)).template get_val<T>();
    }

    const std::string& help() const { return this->help_text; }
    const std::string& get_program_name() const { return this->program_name; }
    const std::string& get_positional(const size_t &pos) const { return this->positionals.at(pos); }
    const bool& is_initialized(const std::string &l_name) const { return this->longargs.at("--" + l_name)->is_initialized(); }
    const bool& is_initialized(const char &l_name) const { return this->shortargs.at(l_name)->is_initialized(); }
    template <typename T> void set_value_ext(const std::string &name, T in_val) {
      this->longargs.at("--" + name)->set_val<T>(in_val);
    }
  };
  static const std::string get_version() { return 'v' + std::to_string(CXXARGS_VERSION_MAJOR) + '.' + std::to_string(CXXARGS_VERSION_MINOR) + '.' + std::to_string(CXXARGS_VERSION_PATCH); }
}

#endif
