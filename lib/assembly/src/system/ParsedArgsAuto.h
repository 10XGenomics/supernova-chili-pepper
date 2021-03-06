///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//File automatically generated by makeParsedArgsAuto.pl
//DO NOT EDIT BY HAND!

#define CommandArgument_StringSet(NAME) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", "", "<required>", "", "");

#define CommandArgument_StringSet_Doc(NAME, DOC) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", "", "<required>", "", DOC);

#define CommandArgument_StringSet_OrDefault(NAME, DEFAULT) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", "", DEFAULT, "", "");

#define CommandArgument_StringSet_OrDefault_Doc(NAME, DEFAULT, DOC) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", "", DEFAULT, "", DOC);

#define CommandArgument_StringSet_Abbr(NAME, ABBR) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", #ABBR, "<required>", "", "");

#define CommandArgument_StringSet_Abbr_Doc(NAME, ABBR, DOC) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", #ABBR, "<required>", "", DOC);

#define CommandArgument_StringSet_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", #ABBR, DEFAULT, "", "");

#define CommandArgument_StringSet_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    vec<String> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<String>", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_Bool(NAME) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", "", "<required>", "", "");

#define CommandArgument_Bool_Doc(NAME, DOC) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", "", "<required>", "", DOC);

#define CommandArgument_Bool_OrDefault(NAME, DEFAULT) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", "", DEFAULT, "", "");

#define CommandArgument_Bool_OrDefault_Doc(NAME, DEFAULT, DOC) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", "", DEFAULT, "", DOC);

#define CommandArgument_Bool_Abbr(NAME, ABBR) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", #ABBR, "<required>", "", "");

#define CommandArgument_Bool_Abbr_Doc(NAME, ABBR, DOC) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", #ABBR, "<required>", "", DOC);

#define CommandArgument_Bool_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", #ABBR, DEFAULT, "", "");

#define CommandArgument_Bool_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    Bool NAME=False; \
    command.ProcessArgs(NAME, #NAME, "Bool", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_LongLongSet(NAME) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", "", "<required>", "", "");

#define CommandArgument_LongLongSet_Doc(NAME, DOC) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", "", "<required>", "", DOC);

#define CommandArgument_LongLongSet_OrDefault(NAME, DEFAULT) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", "", DEFAULT, "", "");

#define CommandArgument_LongLongSet_OrDefault_Doc(NAME, DEFAULT, DOC) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", "", DEFAULT, "", DOC);

#define CommandArgument_LongLongSet_Abbr(NAME, ABBR) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", #ABBR, "<required>", "", "");

#define CommandArgument_LongLongSet_Abbr_Doc(NAME, ABBR, DOC) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", #ABBR, "<required>", "", DOC);

#define CommandArgument_LongLongSet_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", #ABBR, DEFAULT, "", "");

#define CommandArgument_LongLongSet_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    vec<longlong> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<longlong>", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_Int(NAME) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", "<required>", "", "");

#define CommandArgument_Int_Doc(NAME, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", "<required>", "", DOC);

#define CommandArgument_Int_Valid(NAME, VALID) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", "<required>", VALID, "");

#define CommandArgument_Int_Valid_Doc(NAME, VALID, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", "<required>", VALID, DOC);

#define CommandArgument_Int_OrDefault(NAME, DEFAULT) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", DEFAULT, "", "");

#define CommandArgument_Int_OrDefault_Doc(NAME, DEFAULT, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", DEFAULT, "", DOC);

#define CommandArgument_Int_OrDefault_Valid(NAME, DEFAULT, VALID) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", DEFAULT, VALID, "");

#define CommandArgument_Int_OrDefault_Valid_Doc(NAME, DEFAULT, VALID, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", "", DEFAULT, VALID, DOC);

#define CommandArgument_Int_Abbr(NAME, ABBR) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, "<required>", "", "");

#define CommandArgument_Int_Abbr_Doc(NAME, ABBR, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, "<required>", "", DOC);

#define CommandArgument_Int_Abbr_Valid(NAME, ABBR, VALID) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, "<required>", VALID, "");

#define CommandArgument_Int_Abbr_Valid_Doc(NAME, ABBR, VALID, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, "<required>", VALID, DOC);

#define CommandArgument_Int_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, DEFAULT, "", "");

#define CommandArgument_Int_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_Int_Abbr_OrDefault_Valid(NAME, ABBR, DEFAULT, VALID) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, DEFAULT, VALID, "");

#define CommandArgument_Int_Abbr_OrDefault_Valid_Doc(NAME, ABBR, DEFAULT, VALID, DOC) \
    int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "int", #ABBR, DEFAULT, VALID, DOC);

#define CommandArgument_DoubleSet(NAME) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", "", "<required>", "", "");

#define CommandArgument_DoubleSet_Doc(NAME, DOC) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", "", "<required>", "", DOC);

#define CommandArgument_DoubleSet_OrDefault(NAME, DEFAULT) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", "", DEFAULT, "", "");

#define CommandArgument_DoubleSet_OrDefault_Doc(NAME, DEFAULT, DOC) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", "", DEFAULT, "", DOC);

#define CommandArgument_DoubleSet_Abbr(NAME, ABBR) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", #ABBR, "<required>", "", "");

#define CommandArgument_DoubleSet_Abbr_Doc(NAME, ABBR, DOC) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", #ABBR, "<required>", "", DOC);

#define CommandArgument_DoubleSet_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", #ABBR, DEFAULT, "", "");

#define CommandArgument_DoubleSet_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    vec<double> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<double>", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_String(NAME) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", "<required>", "", "");

#define CommandArgument_String_Doc(NAME, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", "<required>", "", DOC);

#define CommandArgument_String_Valid(NAME, VALID) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", "<required>", VALID, "");

#define CommandArgument_String_Valid_Doc(NAME, VALID, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", "<required>", VALID, DOC);

#define CommandArgument_String_OrDefault(NAME, DEFAULT) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", DEFAULT, "", "");

#define CommandArgument_String_OrDefault_Doc(NAME, DEFAULT, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", DEFAULT, "", DOC);

#define CommandArgument_String_OrDefault_Valid(NAME, DEFAULT, VALID) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", DEFAULT, VALID, "");

#define CommandArgument_String_OrDefault_Valid_Doc(NAME, DEFAULT, VALID, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", "", DEFAULT, VALID, DOC);

#define CommandArgument_String_Abbr(NAME, ABBR) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, "<required>", "", "");

#define CommandArgument_String_Abbr_Doc(NAME, ABBR, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, "<required>", "", DOC);

#define CommandArgument_String_Abbr_Valid(NAME, ABBR, VALID) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, "<required>", VALID, "");

#define CommandArgument_String_Abbr_Valid_Doc(NAME, ABBR, VALID, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, "<required>", VALID, DOC);

#define CommandArgument_String_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, DEFAULT, "", "");

#define CommandArgument_String_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_String_Abbr_OrDefault_Valid(NAME, ABBR, DEFAULT, VALID) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, DEFAULT, VALID, "");

#define CommandArgument_String_Abbr_OrDefault_Valid_Doc(NAME, ABBR, DEFAULT, VALID, DOC) \
    String NAME; \
    command.ProcessArgs(NAME, #NAME, "String", #ABBR, DEFAULT, VALID, DOC);

#define CommandArgument_LongLong(NAME) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", "<required>", "", "");

#define CommandArgument_LongLong_Doc(NAME, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", "<required>", "", DOC);

#define CommandArgument_LongLong_Valid(NAME, VALID) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", "<required>", VALID, "");

#define CommandArgument_LongLong_Valid_Doc(NAME, VALID, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", "<required>", VALID, DOC);

#define CommandArgument_LongLong_OrDefault(NAME, DEFAULT) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", DEFAULT, "", "");

#define CommandArgument_LongLong_OrDefault_Doc(NAME, DEFAULT, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", DEFAULT, "", DOC);

#define CommandArgument_LongLong_OrDefault_Valid(NAME, DEFAULT, VALID) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", DEFAULT, VALID, "");

#define CommandArgument_LongLong_OrDefault_Valid_Doc(NAME, DEFAULT, VALID, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", "", DEFAULT, VALID, DOC);

#define CommandArgument_LongLong_Abbr(NAME, ABBR) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, "<required>", "", "");

#define CommandArgument_LongLong_Abbr_Doc(NAME, ABBR, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, "<required>", "", DOC);

#define CommandArgument_LongLong_Abbr_Valid(NAME, ABBR, VALID) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, "<required>", VALID, "");

#define CommandArgument_LongLong_Abbr_Valid_Doc(NAME, ABBR, VALID, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, "<required>", VALID, DOC);

#define CommandArgument_LongLong_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, DEFAULT, "", "");

#define CommandArgument_LongLong_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_LongLong_Abbr_OrDefault_Valid(NAME, ABBR, DEFAULT, VALID) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, DEFAULT, VALID, "");

#define CommandArgument_LongLong_Abbr_OrDefault_Valid_Doc(NAME, ABBR, DEFAULT, VALID, DOC) \
    longlong NAME=0; \
    command.ProcessArgs(NAME, #NAME, "longlong", #ABBR, DEFAULT, VALID, DOC);

#define CommandArgument_UnsignedInt(NAME) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", "<required>", "", "");

#define CommandArgument_UnsignedInt_Doc(NAME, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", "<required>", "", DOC);

#define CommandArgument_UnsignedInt_Valid(NAME, VALID) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", "<required>", VALID, "");

#define CommandArgument_UnsignedInt_Valid_Doc(NAME, VALID, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", "<required>", VALID, DOC);

#define CommandArgument_UnsignedInt_OrDefault(NAME, DEFAULT) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", DEFAULT, "", "");

#define CommandArgument_UnsignedInt_OrDefault_Doc(NAME, DEFAULT, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", DEFAULT, "", DOC);

#define CommandArgument_UnsignedInt_OrDefault_Valid(NAME, DEFAULT, VALID) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", DEFAULT, VALID, "");

#define CommandArgument_UnsignedInt_OrDefault_Valid_Doc(NAME, DEFAULT, VALID, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", "", DEFAULT, VALID, DOC);

#define CommandArgument_UnsignedInt_Abbr(NAME, ABBR) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, "<required>", "", "");

#define CommandArgument_UnsignedInt_Abbr_Doc(NAME, ABBR, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, "<required>", "", DOC);

#define CommandArgument_UnsignedInt_Abbr_Valid(NAME, ABBR, VALID) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, "<required>", VALID, "");

#define CommandArgument_UnsignedInt_Abbr_Valid_Doc(NAME, ABBR, VALID, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, "<required>", VALID, DOC);

#define CommandArgument_UnsignedInt_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, DEFAULT, "", "");

#define CommandArgument_UnsignedInt_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_UnsignedInt_Abbr_OrDefault_Valid(NAME, ABBR, DEFAULT, VALID) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, DEFAULT, VALID, "");

#define CommandArgument_UnsignedInt_Abbr_OrDefault_Valid_Doc(NAME, ABBR, DEFAULT, VALID, DOC) \
    unsigned int NAME=0; \
    command.ProcessArgs(NAME, #NAME, "unsigned int", #ABBR, DEFAULT, VALID, DOC);

#define CommandArgument_Double(NAME) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", "<required>", "", "");

#define CommandArgument_Double_Doc(NAME, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", "<required>", "", DOC);

#define CommandArgument_Double_Valid(NAME, VALID) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", "<required>", VALID, "");

#define CommandArgument_Double_Valid_Doc(NAME, VALID, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", "<required>", VALID, DOC);

#define CommandArgument_Double_OrDefault(NAME, DEFAULT) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", DEFAULT, "", "");

#define CommandArgument_Double_OrDefault_Doc(NAME, DEFAULT, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", DEFAULT, "", DOC);

#define CommandArgument_Double_OrDefault_Valid(NAME, DEFAULT, VALID) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", DEFAULT, VALID, "");

#define CommandArgument_Double_OrDefault_Valid_Doc(NAME, DEFAULT, VALID, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", "", DEFAULT, VALID, DOC);

#define CommandArgument_Double_Abbr(NAME, ABBR) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, "<required>", "", "");

#define CommandArgument_Double_Abbr_Doc(NAME, ABBR, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, "<required>", "", DOC);

#define CommandArgument_Double_Abbr_Valid(NAME, ABBR, VALID) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, "<required>", VALID, "");

#define CommandArgument_Double_Abbr_Valid_Doc(NAME, ABBR, VALID, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, "<required>", VALID, DOC);

#define CommandArgument_Double_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, DEFAULT, "", "");

#define CommandArgument_Double_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, DEFAULT, "", DOC);

#define CommandArgument_Double_Abbr_OrDefault_Valid(NAME, ABBR, DEFAULT, VALID) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, DEFAULT, VALID, "");

#define CommandArgument_Double_Abbr_OrDefault_Valid_Doc(NAME, ABBR, DEFAULT, VALID, DOC) \
    double NAME=0; \
    command.ProcessArgs(NAME, #NAME, "double", #ABBR, DEFAULT, VALID, DOC);

#define CommandArgument_IntSet(NAME) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", "", "<required>", "", "");

#define CommandArgument_IntSet_Doc(NAME, DOC) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", "", "<required>", "", DOC);

#define CommandArgument_IntSet_OrDefault(NAME, DEFAULT) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", "", DEFAULT, "", "");

#define CommandArgument_IntSet_OrDefault_Doc(NAME, DEFAULT, DOC) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", "", DEFAULT, "", DOC);

#define CommandArgument_IntSet_Abbr(NAME, ABBR) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", #ABBR, "<required>", "", "");

#define CommandArgument_IntSet_Abbr_Doc(NAME, ABBR, DOC) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", #ABBR, "<required>", "", DOC);

#define CommandArgument_IntSet_Abbr_OrDefault(NAME, ABBR, DEFAULT) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", #ABBR, DEFAULT, "", "");

#define CommandArgument_IntSet_Abbr_OrDefault_Doc(NAME, ABBR, DEFAULT, DOC) \
    vec<int> NAME; \
    command.ProcessArgs(NAME, #NAME, "vec<int>", #ABBR, DEFAULT, "", DOC);

