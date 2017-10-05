gemma = function(arguments) {
	args = as.character(arguments)
	.C("Rmain",as.integer(length(args)),as.integer(max(nchar(args))),args,EXIT_STATUS=integer(1))$EXIT_STATUS
}
