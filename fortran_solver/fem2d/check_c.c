#include <iostream>

// Declare the external subroutine
extern "C" {
    //void generateellipse_();
    void sayhello();
}

int main() {
    // Call the subroutine
    //generateellipse_();
    sayhello();

    return 0;
}
