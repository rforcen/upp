// the NC zlib compression facilities

#pragma once

#include "numCpp.h"
#include <zlib.h>

template <typename T>
class NCzlib {
 public:
  NCzlib(NC<T> &a) {
    vc = compressVector(a.data);  // compress & release a.data
    a.data.clear();
    ta = &a;
  }
  void decompress() {
    ta->data = decompressVector(vc);  // decompress & release vc
    vc.clear();
  }

  int size() { return vc.size(); }

 private:
  vector<Byte> vc;      // stores compressed data
  NC<T> *ta = nullptr;  // temp NC

 private:
  vector<Byte> compressVector(const vector<T> &data) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, Z_BEST_COMPRESSION) != Z_OK) {
      throw(runtime_error("deflateInit failed while compressing."));
    }

    zs.next_in = reinterpret_cast<Byte *>(const_cast<T *>(data.data()));
    zs.avail_in = data.size() * sizeof(T);

    int ret;
    char outbuffer[32768];
    vector<Byte> outdata;

    do {
      zs.next_out = reinterpret_cast<Byte *>(outbuffer);
      zs.avail_out = sizeof(outbuffer);

      ret = deflate(&zs, Z_FINISH);

      if (outdata.size() < zs.total_out) {
        outdata.insert(outdata.end(), outbuffer,
                       outbuffer + zs.total_out - outdata.size());
      }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {
      throw(runtime_error("Exception during zlib compression: " +
                          to_string(ret)));
    }

    return outdata;
  }

  vector<T> decompressVector(const vector<Byte> &compressedData) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK) 
      throw(runtime_error("d1. inflateInit failed while decompressing."));
    

    zs.next_in =
        reinterpret_cast<Byte *>(const_cast<Byte *>(compressedData.data()));
    zs.avail_in = compressedData.size();

    int ret;
    char outbuffer[32768];
    vector<T> outdata;

    do {
      zs.next_out = reinterpret_cast<Byte *>(outbuffer);
      zs.avail_out = sizeof(outbuffer);

      ret = inflate(&zs, 0);

      if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
        inflateEnd(&zs);
        throw(runtime_error("d2. Exception during zlib decompression: " +
                            to_string(ret)));
      }

      size_t bytesDecompressed = sizeof(outbuffer) - zs.avail_out;
      outdata.insert(
          outdata.end(), reinterpret_cast<T *>(outbuffer),
          reinterpret_cast<T *>(outbuffer) + bytesDecompressed / sizeof(T));
    } while (ret != Z_STREAM_END);

    inflateEnd(&zs);

    return outdata;
  }
};