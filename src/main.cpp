#include <iostream>

#include "foldalign.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Foldalign fold_instance("CCCCCAAAAGGGGG", "CCAACAUUAGUUGG");

    return fold_instance.fold_align();
}
