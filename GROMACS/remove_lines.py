import argparse

def remove_lines_starting_with(file_path, prefixes = ('@', '#')):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        filtered_lines = [line for line in lines if not line.startswith(tuple(prefixes))]

        with open(file_path, 'w') as file:
            file.writelines(filtered_lines)

        print(f"Lines starting with {prefixes} removed from {file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

def remove_lines(file_path, num_lines):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        with open(file_path, 'w') as file:
            file.writelines(lines[num_lines:])
        print(f"{num_lines} lines removed from {file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lines', help='lines to remove', required=False, type=int)
    parser.add_argument('--file', help='file to remove lines', required=True)
    args = parser.parse_args()

    file_path = args.file  # Replace with the actual path to your file
    print(file_path)
    prefixes_to_remove = ('@', '#')

    remove_lines_starting_with(file_path, prefixes_to_remove)

    if args.lines:
        num_lines_to_remove = args.lines
        remove_lines(file_path, num_lines_to_remove)

if __name__ == "__main__":
    main()
