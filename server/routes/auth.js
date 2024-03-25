const express = require('express');
const router = express.Router();
const bcrypt = require('bcryptjs');
const jwt = require('jsonwebtoken');
const User = require('../models/user');
const config = require('../config/environment');

// Login Route
router.post('/login', async (req, res) => {
    const { email, password } = req.body;
    try {
      const user = await User.findOne({ email });
      if (!user) {
        return res.status(401).json({ message: 'Authentication failed. User not found.' });
      }
      bcrypt.compare(password, user.password, (err, isMatch) => {
        if (isMatch && !err) {

          tokenExpirationSeconds = 10 *60 // 10 minutes

          const token = jwt.sign({ _id: user._id }, config.jwtSecret, { expiresIn: tokenExpirationSeconds });
          res.cookie('jwt', token, { httpOnly: true, secure: true, sameSite: 'Strict', maxAge: tokenExpirationSeconds * 1000}); 

          
          const userData = { alias: `${user.name} ${user.lastname[0]}.`, _id: user._id.toString() };
          res.cookie('userData', JSON.stringify(userData), { httpOnly: false, secure: true, sameSite: 'Strict', maxAge: tokenExpirationSeconds * 1000 });

          res.json({ message: "User logged in successfully"});
        } else {
          res.status(401).json({ message: 'Authentication failed. Wrong password.' });
        }
      });
    } catch (error) {
      res.status(500).json({ message: 'Internal server error' });
    }
});
  
// Register Route
router.post('/register', async (req, res) => {
    const { email, password, name, lastname } = req.body;
    if (!email || !password || !name || !lastname) {
      return res.status(400).json({ message: 'All fields are required' });
    }
    try {
      const passwordHash = await bcrypt.hash(password, 10);
      const user = await User.create({ email, password:passwordHash, name, lastname });
      res.json({ message: 'User created successfully', _id:user._id });
    } catch (error) {
      res.status(500).json({ message: 'Internal server error' });
    }
});

module.exports = router;
