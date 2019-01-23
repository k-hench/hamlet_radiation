body() {
	IFS= read -r header
	printf '%s\n' "$header"
	"$@"
}