#include "tsp.h"

/*!
 * This function will parse the command line and saves it into instance
 * @param	argc is the number if string pointed by argv
 * @param	argv is a pointer to a pointer of strings
 * @param	inst pointer to an instance
 * @modify	saves the value passed in argv into inst
 */
void parse_command_line(int argc, char **argv, instance *inst);

/*!
* This fuction will read the file saved into inst and saves
* its value into the instance itself
* @param	inst the instance where the name of file is saved and
*					where the the value will be saved
* @modify	saves the value read into the file into inst
*/
void read_input(instance *inst);

