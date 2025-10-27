validate_input = function(input){
    valid_inputs = c("Q","1","2","3","4","5")
    return (toupper(trimws(input)) %in% valid_inputs)
}

clear_console = function(){
    system("clear")
}

is_valid_input_boolean = function(input){
    return(toupper(input) == 'TRUE' || toupper(input) == 'FALSE')
}