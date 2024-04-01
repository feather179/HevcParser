#include <iostream>
#include <memory>

#include "HevcDecoder.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s INPUT_FILE OUTPUT_YUV\n", argv[0]);
        return 0;
    }

    auto decoder = std::make_shared<HevcDecoder>();
    decoder->decode(argv[1], argv[2]);

    return 0;
}
