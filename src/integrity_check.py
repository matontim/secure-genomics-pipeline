import hashlib
import sys

def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python src/integrity_check.py path/to/file")
    else:
        print(sha256_of_file(sys.argv[1]))
